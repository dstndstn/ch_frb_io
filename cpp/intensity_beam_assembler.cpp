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


intensity_beam_assembler::intensity_beam_assembler(const intensity_network_stream::initializer &x, int assembler_ix) :
    _initializer(x)
{
    if ((assembler_ix < 0) || (assembler_ix >= (int)x.beam_ids.size()))
	throw runtime_error("ch_frb_io: bad assembler_ix passed to intensity_beam_assembler constructor");

    this->beam_id = x.beam_ids[assembler_ix];

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_initflag_set, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);
}


intensity_beam_assembler::~intensity_beam_assembler()
{
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_cond_destroy(&this->cond_initflag_set);
    pthread_mutex_destroy(&this->lock);
}


bool intensity_beam_assembler::wait_for_first_packet(int &nupfreq_, int &nt_per_packet_, int &fpga_counts_per_sample_)
{
    pthread_mutex_lock(&this->lock);

    while (!this->initflag_protected)
	pthread_cond_wait(&this->cond_initflag_set, &this->lock);

    bool alive = !doneflag_protected;
    pthread_mutex_unlock(&this->lock);

    nupfreq_ = this->nupfreq;
    nt_per_packet_ = this->nt_per_packet;
    fpga_counts_per_sample_ = this->fpga_counts_per_sample;

    return alive;
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

	if (this->doneflag_protected) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_assembled_chunks_added, &this->lock);
    }
}


// -------------------------------------------------------------------------------------------------
//
// Routines called by network thread


void intensity_beam_assembler::_put_unassembled_packet(const intensity_packet &packet)
{
    uint64_t packet_it0 = packet.fpga_count / packet.fpga_counts_per_sample;
    uint64_t packet_it1 = packet_it0 + packet.ntsamp;

    if (doneflag_unprotected)
	throw runtime_error("ch_frb_io: internal error: intensity_beam_assembler::_put_unassembled_packet() called after _end_stream()");

    if (!initflag_unprotected) {
	uint64_t assembler_it0 = (packet_it0 / constants::nt_per_assembled_chunk) * constants::nt_per_assembled_chunk;

	this->nupfreq = packet.nupfreq;
	this->nt_per_packet = packet.ntsamp;
	this->fpga_counts_per_sample = packet.fpga_counts_per_sample;

	this->active_chunk0 = this->_make_assembled_chunk(assembler_it0);
	this->active_chunk1 = this->_make_assembled_chunk(active_chunk0->chunk_t1);
	this->initflag_unprotected = true;

	pthread_mutex_lock(&this->lock);
	this->initflag_protected = true;
	pthread_cond_broadcast(&this->cond_initflag_set);
	pthread_mutex_unlock(&this->lock);
    }

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

    // FIXME bookkeep drops!

    if ((packet_it0 >= active_chunk0->chunk_t0) && (packet_it1 <= active_chunk0->chunk_t1))
	active_chunk0->add_packet(packet);
    else if ((packet_it0 >= active_chunk1->chunk_t0) && (packet_it1 <= active_chunk1->chunk_t1))
	active_chunk1->add_packet(packet);
    else if ((packet_it0 < active_chunk1->chunk_t1) && (packet_it1 > active_chunk0->chunk_t0))
	throw runtime_error("DOH");
}


void intensity_beam_assembler::_put_assembled_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    pthread_mutex_lock(&this->lock);

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


void intensity_beam_assembler::_end_stream()
{
    if (this->doneflag_unprotected)
	throw runtime_error("ch_frb_io: internal error: double call to intensity_beam_assembler::_end_stream()");

    if (this->initflag_unprotected) {
	this->_put_assembled_chunk(active_chunk0);
	this->_put_assembled_chunk(active_chunk1);
	this->active_chunk0 = this->active_chunk1 = shared_ptr<assembled_chunk> ();
    }

    this->initflag_unprotected = true;
    this->doneflag_unprotected = true;

    pthread_mutex_lock(&this->lock);
    this->initflag_protected = true;
    this->doneflag_protected = true;
    pthread_cond_broadcast(&this->cond_initflag_set);
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);
}


std::shared_ptr<assembled_chunk> intensity_beam_assembler::_make_assembled_chunk(uint64_t chunk_t0)
{
    if (_initializer.mandate_fast_kernels)
	return make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
    else if (_initializer.mandate_reference_kernels)
	return make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
    else
	return assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
}


}  // namespace ch_frb_io
