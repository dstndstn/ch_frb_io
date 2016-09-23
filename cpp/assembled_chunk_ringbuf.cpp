#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


assembled_chunk_ringbuf::assembled_chunk_ringbuf(const intensity_network_stream &s, int assembler_ix) :
    _initializer(s._initializer)
{
    if ((assembler_ix < 0) || (assembler_ix >= s.nassemblers))
	throw runtime_error("ch_frb_io: bad assembler_ix passed to assembled_chunk_ringbuf constructor");

    this->beam_id = s._initializer.beam_ids[assembler_ix];
    this->nupfreq = s.fp_nupfreq;
    this->nt_per_packet = s.fp_nt_per_packet;
    this->fpga_counts_per_sample = s.fp_fpga_counts_per_sample;

    uint64_t packet_it0 = s.fp_fpga_count / fpga_counts_per_sample;
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
bool assembled_chunk_ringbuf::put_unassembled_packet(const intensity_packet &packet)
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
	this->_put_assembled_chunk(active_chunk0);
	active_chunk0 = active_chunk1;
	active_chunk1 = this->_make_assembled_chunk(active_chunk1->chunk_t1);
    }

    if ((packet_it0 >= active_chunk0->chunk_t0) && (packet_it1 <= active_chunk0->chunk_t1))
	active_chunk0->add_packet(packet);
    else if ((packet_it0 >= active_chunk1->chunk_t0) && (packet_it1 <= active_chunk1->chunk_t1))
	active_chunk1->add_packet(packet);
    else if ((packet_it0 < active_chunk1->chunk_t1) && (packet_it1 > active_chunk0->chunk_t0))
	throw runtime_error("DOH");
    else
	return false;

    return true;
}


void assembled_chunk_ringbuf::_put_assembled_chunk(const shared_ptr<assembled_chunk> &chunk)
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
	cerr << "ch_frb_io: warning: assembler's \"downstream\" thread is running too slow, dropping assembled_chunk\n";
	
	if (!_initializer.drops_allowed)
	    throw runtime_error("ch_frb_io: assembled_chunk was dropped and assembler's 'drops_allowed' flag was set to false");

	return;
    }

    int i = (assembled_ringbuf_pos + assembled_ringbuf_size) % constants::assembled_ringbuf_capacity;
    this->assembled_ringbuf[i] = chunk;
    this->assembled_ringbuf_size++;

    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);
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


void assembled_chunk_ringbuf::end_stream()
{
    if (!active_chunk0 || !active_chunk1)
	throw runtime_error("ch_frb_io: internal error: empty pointers in assembled_chunk_ringbuf::end_stream(), this can happen if end_stream() is called twice");

    this->_put_assembled_chunk(active_chunk0);
    this->_put_assembled_chunk(active_chunk1);
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
    if (_initializer.mandate_fast_kernels)
	return make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
    else if (_initializer.mandate_reference_kernels)
	return make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
    else
	return assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
}


}  // namespace ch_frb_io
