#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// class intensity_beam_assembler


// static member function used as de facto constructor
shared_ptr<intensity_beam_assembler> intensity_beam_assembler::make(int beam_id_, bool drops_allowed_)
{
    intensity_beam_assembler *ret_bare = new intensity_beam_assembler(beam_id_, drops_allowed_);
    shared_ptr<intensity_beam_assembler> ret(ret_bare);

    int err = pthread_create(&ret->assembler_thread, NULL, intensity_beam_assembler::assembler_pthread_main, &ret);
    if (err)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_beam_assembler constructor: ") + strerror(errno));

    pthread_mutex_lock(&ret->lock);
    while (!ret->assembler_thread_started)
	pthread_cond_wait(&ret->cond_assembler_state_changed, &ret->lock);
    pthread_mutex_unlock(&ret->lock);

    return ret;
}


intensity_beam_assembler::intensity_beam_assembler(int beam_id_, bool drops_allowed_) 
    : beam_id(beam_id_), drops_allowed(drops_allowed_)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("intensity_beam_constructor: invalid beam_id");

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembler_state_changed, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);

    int capacity = constants::unassembled_ringbuf_capacity;
    int max_npackets = constants::max_unassembled_packets_per_list;
    int max_nbytes = constants::max_unassembled_nbytes_per_list;
    string dropstr = "warning: assembler thread running too slow, dropping packets";

    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (capacity, max_npackets, max_nbytes, dropstr, drops_allowed);
}


intensity_beam_assembler::~intensity_beam_assembler()
{
    pthread_cond_destroy(&this->cond_assembler_state_changed);
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_mutex_destroy(&this->lock);
}


bool intensity_beam_assembler::wait_for_first_packet(int &fpga_counts_per_sample_, int &nupfreq_)
{
    pthread_mutex_lock(&this->lock);

    while (!first_packet_received)
	pthread_cond_wait(&this->cond_assembler_state_changed, &this->lock);

    fpga_counts_per_sample_ = this->fpga_counts_per_sample;
    nupfreq_ = this->nupfreq;

    pthread_mutex_unlock(&this->lock);

    // FIXME logical choice?  Should we just make this return 'void'?
    return (nupfreq != 0);
}


// advances state: stream_ended -> assembler_ended
void intensity_beam_assembler::assembler_thread_end()
{
    pthread_mutex_lock(&this->lock);

    // We set the whole chain of flags, to avoid confusion.
    this->first_packet_received = true;
    this->assembler_thread_ended = true;

    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);

    // Probably redundant but playing it safe...
    this->unassembled_ringbuf->end_stream();
}


void intensity_beam_assembler::join_assembler_thread()
{
    bool call_join_after_releasing_lock = false;

    pthread_mutex_lock(&this->lock);

    while (!this->assembler_thread_ended)
	pthread_cond_wait(&this->cond_assembler_state_changed, &this->lock);

    if (!this->assembler_thread_joined) {
	this->assembler_thread_joined = true;
	call_join_after_releasing_lock = true;
    }

    pthread_mutex_unlock(&this->lock);
	
    if (call_join_after_releasing_lock)
	pthread_join(this->assembler_thread, NULL);
}


bool intensity_beam_assembler::put_unassembled_packets(udp_packet_list &packet_list)
{
    bool is_blocking = false;
    return this->unassembled_ringbuf->producer_put_packet_list(packet_list, is_blocking);
}


void intensity_beam_assembler::put_assembled_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    pthread_mutex_lock(&this->lock);

    if (assembled_ringbuf_size >= constants::assembled_ringbuf_capacity) {
	pthread_mutex_unlock(&this->lock);
	cerr << "ch_frb_io: warning: assembler's \"downstream\" thread is running too slow, some packets will be dropped\n";
	
	if (!drops_allowed) {
	    cerr << "ch_frb_io: assembler's 'drops_allowed' flag was set to false, this will be treated as an error\n";
	    exit(1);
	}

	return;
    }

    int i = (assembled_ringbuf_pos + assembled_ringbuf_size) % constants::assembled_ringbuf_capacity;
    this->assembled_ringbuf[i] = chunk;
    this->assembled_ringbuf_size++;

    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);
}


bool intensity_beam_assembler::get_assembled_chunk(shared_ptr<assembled_chunk> &chunk)
{
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (assembled_ringbuf_size > 0) {
	    int i = assembled_ringbuf_pos % constants::assembled_ringbuf_capacity;
	    chunk = assembled_ringbuf[i];

	    assembled_ringbuf[i] = shared_ptr<assembled_chunk> ();
	    this->assembled_ringbuf_pos++;
	    this->assembled_ringbuf_size--;	    

	    pthread_mutex_unlock(&this->lock);

	    if (!chunk)
		throw runtime_error("ch_frb_io: internal error: unexpected empty pointer in get_assembled_chunk()");

	    return true;
	}

	if (assembler_thread_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_assembled_chunks_added, &this->lock);
    }
}


// -------------------------------------------------------------------------------------------------
//
// Routines called by network thread


void intensity_beam_assembler::_announce_first_packet(int fpga_counts_per_sample_, int nupfreq_)
{
    pthread_mutex_lock(&this->lock);

    if (!assembler_thread_started) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: intensity_beam_assembler::_announce_first_packet() was called, but no assembler thread");
    }

    if (this->first_packet_received) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: double call to intensity_beam_assembler::_announce_first_packet()");
    }

    this->fpga_counts_per_sample = fpga_counts_per_sample_;
    this->nupfreq = nupfreq_;

    this->first_packet_received = true;
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_mutex_unlock(&this->lock);
}


// advances state: stream_started -> stream_ended
void intensity_beam_assembler::_end_stream()
{
    // FIXME comment on this
    pthread_mutex_lock(&this->lock);
    this->first_packet_received = true;
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_mutex_unlock(&this->lock);

    this->unassembled_ringbuf->end_stream();
}



// -------------------------------------------------------------------------------------------------
//
// assembler thread


// static member function
void *intensity_beam_assembler::assembler_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to assembler_pthread_main()");

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<intensity_beam_assembler> *arg = (shared_ptr<intensity_beam_assembler> *) opaque_arg;
    shared_ptr<intensity_beam_assembler> assembler = *arg;

    if (!assembler)
	throw runtime_error("ch_frb_io: internal error: empty pointer passed to assembler_thread_main()");

    try {
	assembler->assembler_thread_main();
    } catch (...) {
	assembler->assembler_thread_end();
	throw;
    }

    assembler->assembler_thread_end();

    return NULL;
}


void intensity_beam_assembler::assembler_thread_main()
{
    pthread_mutex_lock(&this->lock);

    this->assembler_thread_started = true;
    pthread_cond_broadcast(&this->cond_assembler_state_changed);

    while (!this->first_packet_received)
	pthread_cond_wait(&this->cond_assembler_state_changed, &this->lock);

    pthread_mutex_unlock(&this->lock);

    // Sanity check stream params
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("intensity_beam_assembler: bad value of fpga_counts received from stream");
    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("intensity_beam_assembler: bad value of nupfreq received from stream");

    udp_packet_list unassembled_packet_list(constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list);
    intensity_packet packet;

    // The 'initialized' flag refers to 'assembler_it0', 'chunk0_*', and 'chunk1_*'
    bool initialized = false;
    uint64_t assembler_it0 = 0;

    std::shared_ptr<assembled_chunk> chunk0;
    std::shared_ptr<assembled_chunk> chunk1;

    // Outer loop over unassembled packets

    for (;;) {
	bool alive = unassembled_ringbuf->consumer_get_packet_list(unassembled_packet_list);

	if (!alive) {
	    // FIXME can this be improved?
	    if (!initialized)
		throw runtime_error("ch_frb_io: assembler thread failed to receive any packets from network thread");

	    this->put_assembled_chunk(chunk0);
	    this->put_assembled_chunk(chunk1);
	    return;
	}

	for (int ipacket = 0; ipacket < unassembled_packet_list.curr_npackets; ipacket++) {
	    uint8_t *packet_data = unassembled_packet_list.data_start + unassembled_packet_list.packet_offsets[ipacket];
	    int packet_nbytes = unassembled_packet_list.packet_offsets[ipacket+1] - unassembled_packet_list.packet_offsets[ipacket];	    
	    bool well_formed = packet.read(packet_data, packet_nbytes);

	    // FIXME needs comment
	    if (_unlikely(!well_formed || (packet.nbeams != 1) || (packet.beam_ids[0] != this->beam_id)))
		throw runtime_error("ch_frb_io: internal error in assembler thread");

	    uint64_t packet_it0 = packet.fpga_count / fpga_counts_per_sample;
	    uint64_t packet_it1 = packet_it0 + packet.ntsamp;

	    if (!initialized) {
		// round down to multiple of constants::nt_per_assembled_chunk
		assembler_it0 = (packet_it0 / constants::nt_per_assembled_chunk) * constants::nt_per_assembled_chunk;
		initialized = true;

		chunk0 = make_shared<assembled_chunk> (beam_id, nupfreq, fpga_counts_per_sample, assembler_it0);
		chunk1 = make_shared<assembled_chunk> (beam_id, nupfreq, fpga_counts_per_sample, assembler_it0 + constants::nt_per_assembled_chunk);
	    }

	    if (packet_it1 > assembler_it0 + 2 * constants::nt_per_assembled_chunk) {
		//
		// If we receive a packet whose timestamps extend past the range of our current
		// assembly buffer, then we advance the buffer and send an assembled_chunk to the
		// "downstream" thread.
		//
		// A design decision here: for a packet which is far in the future, we advance the 
		// buffer by one assembled_chunk, rather than using the minimum number of advances
		// needed.  This is to avoid a situation where a single rogue packet timestamped
		// in the far future effectively kills the L1 node.
		//
		this->put_assembled_chunk(chunk0);
		assembler_it0 += constants::nt_per_assembled_chunk;

		chunk0 = chunk1;
		chunk1 = make_shared<assembled_chunk> (beam_id, nupfreq, fpga_counts_per_sample, assembler_it0 + constants::nt_per_assembled_chunk);
	    }

	    if ((packet_it0 >= assembler_it0) && (packet_it1 <= assembler_it0 + constants::nt_per_assembled_chunk)) {
		// packet is a subset of chunk0
		int offset = packet_it0 - assembler_it0;
		packet.decode(chunk0->intensity + offset, chunk0->weights + offset, constants::nt_per_assembled_chunk);
		continue;
	    }

	    if ((packet_it0 >= assembler_it0 + constants::nt_per_assembled_chunk) && (packet_it1 <= assembler_it0 + 2 * constants::nt_per_assembled_chunk)) {
		// packet is a subset of chunk1
		int offset = packet_it0 - assembler_it0 - constants::nt_per_assembled_chunk;
		packet.decode(chunk1->intensity + offset, chunk1->weights + offset, constants::nt_per_assembled_chunk);
		continue;
	    }

	    if ((packet_it1 <= assembler_it0) || (packet_it0 >= assembler_it0 + 2 * constants::nt_per_assembled_chunk))
		continue;

	    throw runtime_error("DOH");
	}
    }
}


}  // namespace ch_frb_io
