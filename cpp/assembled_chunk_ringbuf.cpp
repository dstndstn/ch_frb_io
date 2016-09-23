#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


assembled_chunk_ringbuf::assembled_chunk_ringbuf(const intensity_network_stream::initializer &ini_params_, int beam_id_, int nupfreq_,
						 int nt_per_packet_, uint64_t fpga_counts_per_sample_, uint64_t fpga_count0) :
    ini_params(ini_params_),
    beam_id(beam_id_),
    nupfreq(nupfreq_),
    nt_per_packet(nt_per_packet_),
    fpga_counts_per_sample(fpga_counts_per_sample_)
{
    // FIXME I suppose we should do some parameter checking

    uint64_t packet_it0 = fpga_count0 / fpga_counts_per_sample;
    uint64_t assembler_it0 = (packet_it0 / constants::nt_per_assembled_chunk) * constants::nt_per_assembled_chunk;

    this->active_chunk0 = this->_make_assembled_chunk(assembler_it0);
    this->active_chunk1 = this->_make_assembled_chunk(active_chunk0->chunk_t1);

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);
}


assembled_chunk_ringbuf::~assembled_chunk_ringbuf()
{
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_mutex_destroy(&this->lock);
}


// Returns true if packet successfully assembled
void assembled_chunk_ringbuf::put_unassembled_packet(const intensity_packet &packet, int64_t *event_counts)
{
    if (!active_chunk0 || !active_chunk1)
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");

    uint64_t packet_it0 = packet.fpga_count / packet.fpga_counts_per_sample;
    uint64_t packet_it1 = packet_it0 + packet.ntsamp;

    if (packet_it1 > active_chunk1->chunk_t1) {
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
	this->_put_assembled_chunk(active_chunk0, event_counts);
	active_chunk0 = active_chunk1;
	active_chunk1 = this->_make_assembled_chunk(active_chunk1->chunk_t1);
    }

    if ((packet_it0 >= active_chunk0->chunk_t0) && (packet_it1 <= active_chunk0->chunk_t1)) {
	event_counts[intensity_network_stream::event_type::assembler_hit]++;
	active_chunk0->add_packet(packet);
    }
    else if ((packet_it0 >= active_chunk1->chunk_t0) && (packet_it1 <= active_chunk1->chunk_t1)) {
	event_counts[intensity_network_stream::event_type::assembler_hit]++;
	active_chunk1->add_packet(packet);
    }
    else if ((packet_it1 <= active_chunk0->chunk_t0) || (packet_it0 >= active_chunk1->chunk_t1))
	event_counts[intensity_network_stream::event_type::assembler_miss]++;
    else
	throw runtime_error("ch_frb_io: internal error: bad packet alignment in assembled_chunk_ringbuf::put_unassembled_packet()");
}


void assembled_chunk_ringbuf::_put_assembled_chunk(const shared_ptr<assembled_chunk> &chunk, int64_t *event_counts)
{
    if (!chunk)
	throw runtime_error("ch_frb_io: internal error: empty pointer passed to assembled_chunk_ringbuf::_put_unassembled_packet()");
	
    pthread_mutex_lock(&this->lock);

    if (this->doneflag) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");
    }

    if (assembled_ringbuf_size >= constants::assembled_ringbuf_capacity) {
	pthread_mutex_unlock(&this->lock);
	event_counts[intensity_network_stream::event_type::assembled_chunk_dropped]++;

	cerr << "ch_frb_io: warning: assembler's \"downstream\" thread is running too slow, dropping assembled_chunk\n";
	
	if (ini_params.drops_allowed)
	    throw runtime_error("ch_frb_io: assembled_chunk was dropped and assembler's 'drops_allowed' flag was set to false");

	return;
    }

    int i = (assembled_ringbuf_pos + assembled_ringbuf_size) % constants::assembled_ringbuf_capacity;
    this->assembled_ringbuf[i] = chunk;
    this->assembled_ringbuf_size++;

    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);
    event_counts[intensity_network_stream::event_type::assembled_chunk_queued]++;
}


shared_ptr<assembled_chunk> assembled_chunk_ringbuf::get_assembled_chunk()
{
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (assembled_ringbuf_size > 0) {
	    int i = assembled_ringbuf_pos % constants::assembled_ringbuf_capacity;
	    auto chunk = assembled_ringbuf[i];
	    assembled_ringbuf[i] = shared_ptr<assembled_chunk> ();

	    this->assembled_ringbuf_pos++;
	    this->assembled_ringbuf_size--;
	    pthread_mutex_unlock(&this->lock);

	    if (!chunk)
		throw runtime_error("ch_frb_io: internal error: unexpected empty pointer in get_assembled_chunk()");

	    return chunk;
	}

	if (this->doneflag) {
	    pthread_mutex_unlock(&this->lock);
	    return shared_ptr<assembled_chunk>();
	}

	pthread_cond_wait(&this->cond_assembled_chunks_added, &this->lock);
    }
}


void assembled_chunk_ringbuf::end_stream(int64_t *event_counts)
{
    if (!active_chunk0 || !active_chunk1)
	throw runtime_error("ch_frb_io: internal error: empty pointers in assembled_chunk_ringbuf::end_stream(), this can happen if end_stream() is called twice");

    this->_put_assembled_chunk(active_chunk0, event_counts);
    this->_put_assembled_chunk(active_chunk1, event_counts);
    this->active_chunk0 = this->active_chunk1 = shared_ptr<assembled_chunk> ();

    pthread_mutex_lock(&this->lock);

    if (doneflag) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: doneflag already set in assembled_chunk_ringbuf::end_stream()");
    }

    this->doneflag = true;
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);
}


std::shared_ptr<assembled_chunk> assembled_chunk_ringbuf::_make_assembled_chunk(uint64_t chunk_t0)
{
    if (ini_params.mandate_fast_kernels)
	return make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
    else if (ini_params.mandate_reference_kernels)
	return make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
    else
	return assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
}


}  // namespace ch_frb_io
