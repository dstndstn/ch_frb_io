#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// defined later in this file
static void *assembler_thread_main(void *opaque_arg);
static void  assembler_thread_main2(intensity_beam_assembler *assembler);


// -------------------------------------------------------------------------------------------------
//
// class intensity_beam_assembler


// static member function used as de facto constructor
shared_ptr<intensity_beam_assembler> intensity_beam_assembler::make(int beam_id_)
{
    intensity_beam_assembler *ret_bare = new intensity_beam_assembler(beam_id_);
    shared_ptr<intensity_beam_assembler> ret(ret_bare);

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<intensity_beam_assembler> *arg = new shared_ptr<intensity_beam_assembler> (ret);
    int err = pthread_create(&ret->assembler_thread, NULL, assembler_thread_main, arg);

    if (err) {
        delete arg;   // note: if pthread_create() succeeds, then assembler thread will delete this pointer
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_beam_assembler constructor: ") + strerror(errno));
    }

    ret->wait_for_assembler_thread_startup();
    return ret;
}


intensity_beam_assembler::intensity_beam_assembler(int beam_id_) 
    : beam_id(beam_id_)
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

    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (capacity, max_npackets, max_nbytes, dropstr);
}


intensity_beam_assembler::~intensity_beam_assembler()
{
    pthread_cond_destroy(&this->cond_assembler_state_changed);
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_mutex_destroy(&this->lock);
}


void intensity_beam_assembler::assembler_thread_startup()
{
    pthread_mutex_lock(&this->lock);
    this->assembler_thread_started = true;
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_mutex_unlock(&this->lock);
}


void intensity_beam_assembler::wait_for_assembler_thread_startup()
{
    pthread_mutex_lock(&this->lock);

    while (!assembler_thread_started)
	pthread_cond_wait(&this->cond_assembler_state_changed, &this->lock);

    pthread_mutex_unlock(&this->lock);
}


void intensity_beam_assembler::start_stream(int fpga_counts_per_sample_, int nupfreq_)
{
    pthread_mutex_lock(&this->lock);
    
    if (stream_started) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: double call to intensity_beam_assembler::start_stream()");
    }

    this->fpga_counts_per_sample = fpga_counts_per_sample_;
    this->nupfreq = nupfreq_;
    this->stream_started = true;

    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_mutex_unlock(&this->lock);
}


bool intensity_beam_assembler::wait_for_stream_params(int &fpga_counts_per_sample_, int &nupfreq_)
{
    pthread_mutex_lock(&this->lock);

    while (!stream_started)
	pthread_cond_wait(&this->cond_assembler_state_changed, &this->lock);

    fpga_counts_per_sample_ = this->fpga_counts_per_sample;
    nupfreq_ = this->nupfreq;

    bool retval = !stream_ended;
    pthread_mutex_unlock(&this->lock);
    return retval;
}


void intensity_beam_assembler::end_stream(bool join_thread)
{
    bool call_join_after_releasing_lock = false;

    pthread_mutex_lock(&this->lock);

    // We set both flags to imitate a stream which has been started and stopped.
    // This ensures that threads sleeping in wait_for_stream_params() get unblocked.
    this->stream_started = true;
    this->stream_ended = true;

    // Wake up all threads in sight (overkill but that's OK)
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);

    if (join_thread && !this->assembler_thread_joined) {
	this->assembler_thread_joined = true;
	call_join_after_releasing_lock = true;
    }

    pthread_mutex_unlock(&this->lock);

    this->unassembled_ringbuf->end_stream();

    if (call_join_after_releasing_lock)
	pthread_join(this->assembler_thread, NULL);
}


udp_packet_list intensity_beam_assembler::allocate_unassembled_packet_list() const
{
    return this->unassembled_ringbuf->allocate_packet_list();
}

bool intensity_beam_assembler::put_unassembled_packets(udp_packet_list &packet_list)
{
    bool is_blocking = false;
    return this->unassembled_ringbuf->producer_put_packet_list(packet_list, is_blocking);
}


bool intensity_beam_assembler::get_unassembled_packets(udp_packet_list &packet_list)
{
    return this->unassembled_ringbuf->consumer_get_packet_list(packet_list);
}


void intensity_beam_assembler::put_assembled_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    pthread_mutex_lock(&this->lock);

    if (assembled_ringbuf_size >= constants::assembled_ringbuf_capacity) {
	pthread_mutex_unlock(&this->lock);
	cerr << "ch_frb_io: warning: assembler's \"downstream\" thread is running too slow, some packets will be dropped\n";
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
	    chunk = assembled_ringbuf[assembled_ringbuf_pos % constants::assembled_ringbuf_capacity];
	    this->assembled_ringbuf_pos++;
	    this->assembled_ringbuf_size--;
	    
	    pthread_mutex_unlock(&this->lock);
	    return true;
	}

	if (stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_assembled_chunks_added, &this->lock);
    }
}


// -------------------------------------------------------------------------------------------------
//
// assembler thread


// FIXME planned optimization for full CHIME: write assembly language kernel for special case nt=16
inline void assemble_packet_row(int nupfreq, int nt, float *dst_intensity, float *dst_weights, const uint8_t *src, int src_stride, float scale, float offset)
{
    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	float *dst2_intensity = dst_intensity + iupfreq * constants::nt_per_assembled_chunk;
	float *dst2_weights = dst_weights + iupfreq * constants::nt_per_assembled_chunk;
	const uint8_t *src2 = src + iupfreq * src_stride;

	for (int it = 0; it < nt; it++) {
	    float t = (float)src2[it];
	    dst2_intensity[it] = scale*t + offset;
	    dst2_weights[it] = ((t*(255.-t)) > 0.5) ? 1.0 : 0.0;   // fastest way to compute weight without using assembly language kernel?
	}
    }
}


inline void assemble_packet(int nfreq, const uint16_t *freq_ids, int nupfreq, int nt, 
			    float *dst_intensity, float *dst_weights, const uint8_t *src, 
			    int src_stride, const float *scales, const float *offsets)
{
    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	// Note: freq_id not range-checked here, since already checked in assembler_thread_main().
	int freq_id = (int)freq_ids[ifreq];

	assemble_packet_row(nupfreq, nt, 
			    dst_intensity + freq_id * nupfreq * constants::nt_per_assembled_chunk,
			    dst_weights + freq_id * nupfreq * constants::nt_per_assembled_chunk,
			    src + ifreq * nupfreq * src_stride,
			    src_stride, scales[ifreq], offsets[ifreq]);
    }
}


static void *assembler_thread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to assembler_thread_main()");

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<intensity_beam_assembler> *arg = (shared_ptr<intensity_beam_assembler> *) opaque_arg;
    shared_ptr<intensity_beam_assembler> assembler = *arg;   // 'exchanger' is safe to use below
    delete arg;

    if (!assembler)
	throw runtime_error("ch_frb_io: internal error: empty pointer passed to assembler_thread_main()");

    cerr << ("ch_frb_io: assembler thread starting (beam_id=" + to_string(assembler->beam_id) + ")\n");
    assembler->assembler_thread_startup();

    try {
	assembler_thread_main2(assembler.get());
    } catch (...) {
	assembler->end_stream(false);  // the "false" means "nonblocking" (otherwise we'd deadlock!)
	throw;
    }

    assembler->end_stream(false);  // "false" means same as above
    cerr << ("ch_frb_io: assembler thread ending (beam_id=" + to_string(assembler->beam_id) + ")\n");    

    return NULL;
}


static void assembler_thread_main2(intensity_beam_assembler *assembler)
{
    int assembler_beam_id = assembler->beam_id;
    udp_packet_list unassembled_packet_list = assembler->allocate_unassembled_packet_list();

    // The 'initialized' flag refers to 'assembler_it0', 'chunk0_*', and 'chunk1_*'
    bool initialized = false;
    uint64_t assembler_it0 = 0;

    std::shared_ptr<assembled_chunk> chunk0;
    float *chunk0_intensity = nullptr;
    float *chunk0_weights = nullptr;

    std::shared_ptr<assembled_chunk> chunk1;
    float *chunk1_intensity = nullptr;
    float *chunk1_weights = nullptr;

    int nupfreq = 0;
    int fpga_counts_per_sample = 0;
    if (!assembler->wait_for_stream_params(fpga_counts_per_sample, nupfreq))
	return;

    // Sanity check stream params
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample >= 65536))
	throw runtime_error("intensity_beam_assembler: bad value of fpga_counts received from stream");
    if ((nupfreq <= 0) || (nupfreq > 16))
	throw runtime_error("intensity_beam_assembler: bad value of nupfreq received from stream");

    // Outer loop over unassembled packets

    for (;;) {
	bool alive = assembler->get_unassembled_packets(unassembled_packet_list);

	if (!alive) {
	    // FIXME can this be improved?
	    if (!initialized)
		throw runtime_error("ch_frb_io: assembler thread failed to receive any packets from network thread");

	    cerr << ("ch_frb_io: assembler thread input is complete (beam_id=" + to_string(assembler_beam_id) + ")\n");
	    assembler->put_assembled_chunk(chunk0);
	    assembler->put_assembled_chunk(chunk1);
	    return;
	}

	for (int ipacket = 0; ipacket < unassembled_packet_list.curr_npackets; ipacket++) {
	    const uint8_t *packet = unassembled_packet_list.data_start + unassembled_packet_list.packet_offsets[ipacket];
	    int packet_nbytes = unassembled_packet_list.packet_offsets[ipacket+1] - unassembled_packet_list.packet_offsets[ipacket];

	    if (_unlikely(packet_nbytes < 24))
		continue;  // FIXME skip

	    int data_nbytes = *((int16_t *) (packet+4));
	    int packet_fpga_counts_per_sample = *((uint16_t *) (packet+6));
	    uint64_t packet_fpga_count = *((uint64_t *) (packet+8));
	    int packet_nbeam = *((uint16_t *) (packet+16));
	    int packet_nfreq = *((uint16_t *) (packet+18));
	    int packet_nupfreq = *((uint16_t *) (packet+20));
	    int packet_ntsamp = *((uint16_t *) (packet+22));

	    bool bad_packet = false;
	    bad_packet |= (packet_fpga_counts_per_sample != fpga_counts_per_sample);
	    bad_packet |= (packet_fpga_count % fpga_counts_per_sample);
	    bad_packet |= (packet_nbeam != 1);
	    bad_packet |= (packet_nfreq <= 0);
	    bad_packet |= (packet_nupfreq != nupfreq);
	    bad_packet |= (packet_ntsamp <= 0);
	    bad_packet |= (packet_ntsamp > constants::nt_per_assembled_chunk);
	    bad_packet |= (data_nbytes != (packet_nfreq * packet_nupfreq * packet_ntsamp));
	    bad_packet |= (packet_nbytes != packet_size(1,packet_nfreq,nupfreq,packet_ntsamp));

	    if (_unlikely(bad_packet))
		continue;  // FIXME skip

	    // Note: all time indices in this function (and in the assembled_chunk struct) are
	    // in units of "downsampled intensities", i.e. (fpga counts) / fpga_counts_per_sample.
	    uint64_t packet_it0 = packet_fpga_count / fpga_counts_per_sample;

	    const uint16_t *packet_beam_ids = (const uint16_t *) (packet + 24);
	    const uint16_t *packet_freq_ids = (const uint16_t *) (packet + 26);
	    const float *packet_scales = (const float *) (packet + 26 + 2*packet_nfreq);
	    const float *packet_offsets = (const float *) (packet + 26 + 6*packet_nfreq);
	    const uint8_t *packet_data = (packet + 26 + 10*packet_nfreq);

	    bad_packet = false;
	    bad_packet |= (packet_beam_ids[0] != assembler_beam_id);
	    
	    for (int ifreq = 0; ifreq < packet_nfreq; ifreq++) {
		int freq_id = packet_freq_ids[ifreq];
		bad_packet |= (freq_id < 0);
		bad_packet |= (freq_id >= 1024);   // FIXME hardcoded constant here
	    }

	    if (_unlikely(bad_packet))
		continue;  // FIXME skip

	    if (!initialized) {
		// round down to multiple of constants::nt_per_assembled_chunk
		assembler_it0 = (packet_it0 / constants::nt_per_assembled_chunk) * constants::nt_per_assembled_chunk;
		initialized = true;

		chunk0 = make_shared<assembled_chunk> (assembler_beam_id, nupfreq, fpga_counts_per_sample, assembler_it0);
		chunk0_intensity = chunk0->intensity;
		chunk0_weights = chunk0->weights;

		chunk1 = make_shared<assembled_chunk> (assembler_beam_id, nupfreq, fpga_counts_per_sample, assembler_it0 + constants::nt_per_assembled_chunk);
		chunk1_intensity = chunk1->intensity;
		chunk1_weights = chunk1->weights;
	    }

	    if (packet_it0 + packet_ntsamp > assembler_it0 + 2 * constants::nt_per_assembled_chunk) {
		//
		// If we receive a packet whose timestamps extend past the range of our current
		// assembly buffer, then we advance the buffer and send an assembled_chunk to the
		// "downstream" thread.
		//
		// A design decision here: for a packet which is far in the future, we advance the 
		// buffer by one assembled_chunk, rather than using the minimum number of advances
		// needed.  This is to avoid a situation where a single rogue packet timestamped
		// in the far future kills the L1 node.
		//
		assembler->put_assembled_chunk(chunk0);
		assembler_it0 += constants::nt_per_assembled_chunk;

		chunk0 = chunk1;
		chunk0_intensity = chunk1_intensity;
		chunk0_weights = chunk1_weights;

		chunk1 = make_shared<assembled_chunk> (assembler_beam_id, nupfreq, fpga_counts_per_sample, assembler_it0 + constants::nt_per_assembled_chunk);
		chunk1_intensity = chunk1->intensity;
		chunk1_weights = chunk1->weights;		
	    }

	    // Compute packet overlap with first assembled_chunk.
	    uint64_t overlap_it0 = max(packet_it0, assembler_it0);
	    uint64_t overlap_it1 = min(packet_it0 + packet_ntsamp, assembler_it0 + constants::nt_per_assembled_chunk);
	    
	    // If there is an overlap then assemble.
	    if (overlap_it0 < overlap_it1) {
		int dst_dt = overlap_it0 - assembler_it0;
		int src_dt = overlap_it0 - packet_it0;

		assemble_packet(packet_nfreq, packet_freq_ids, nupfreq,
				(int)(overlap_it1 - overlap_it0),   // nt
				chunk0_intensity + dst_dt,          // dst_intensity
				chunk0_weights + dst_dt,            // dst_weights
				packet_data + src_dt,               // src
				packet_ntsamp,                      // src_strides
				packet_scales, packet_offsets);
	    }

	    // Logic for the second assembled_chunk follows.
	    overlap_it0 = max(packet_it0, assembler_it0 + constants::nt_per_assembled_chunk);
	    overlap_it1 = min(packet_it0 + packet_ntsamp, assembler_it0 + 2 * constants::nt_per_assembled_chunk);
	    
	    if (overlap_it0 < overlap_it1) {
		int dst_dt = overlap_it0 - (assembler_it0 + constants::nt_per_assembled_chunk);
		int src_dt = overlap_it0 - packet_it0;

		assemble_packet(packet_nfreq, packet_freq_ids, nupfreq,
				(int)(overlap_it1 - overlap_it0),   // nt
				chunk1_intensity + dst_dt,          // dst_intensity
				chunk1_weights + dst_dt,            // dst_weights
				packet_data + src_dt,               // src
				packet_ntsamp,                      // src_strides
				packet_scales, packet_offsets);
	    }	    
	}
    }
}


}  // namespace ch_frb_io
