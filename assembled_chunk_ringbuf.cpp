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
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("ch_frb_io: bad beam_id passed to assembled_chunk_ringbuf constructor");
    if ((nupfreq < 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("ch_frb_io: bad nupfreq value passed to assembled_chunk_ringbuf constructor");
    if ((nt_per_packet <= 0) || (constants::nt_per_assembled_chunk % nt_per_packet != 0))
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::nt_per_packet must be a divisor of constants::nt_per_assembled_chunk");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("ch_frb_io: bad fpga_counts_per_sample value passed to assembled_chunk_ringbuf constructor");
    if (fpga_count0 % fpga_counts_per_sample != 0)
	throw runtime_error("ch_frb_io: assembled_chunk_ringbuf constructor: fpga_count0 was not a multiple of fpga_counts_per_sample");

#ifndef __AVX2__
    if (ini_params.mandate_fast_kernels)
	throw runtime_error("ch_frb_io: the 'mandate_fast_kernels' flag was set, but this machine does not have the AVX2 instruction set");
#endif

    uint64_t packet_t0 = fpga_count0 / fpga_counts_per_sample;
    uint64_t ichunk = packet_t0 / constants::nt_per_assembled_chunk;

    this->active_chunk0 = this->_make_assembled_chunk(ichunk);
    this->active_chunk1 = this->_make_assembled_chunk(ichunk+1);
    this->assembled_ringbuf_pos = ichunk;
    this->assembled_ringbuf_size = 0;

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);
}


assembled_chunk_ringbuf::~assembled_chunk_ringbuf()
{
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_mutex_destroy(&this->lock);
}

int assembled_chunk_ringbuf::get_assembled_ringbuf_size()
{
    int rtn = 0;
    pthread_mutex_lock(&this->lock);
    rtn = this->assembled_ringbuf_size;
    pthread_mutex_unlock(&this->lock);
    return rtn;
}

std::vector<std::shared_ptr<assembled_chunk> >
assembled_chunk_ringbuf::get_ringbuf_snapshot()
{
    std::vector<std::shared_ptr<assembled_chunk> > ring;
    pthread_mutex_lock(&this->lock);
    for (int i=0; i<constants::assembled_ringbuf_capacity; i++) {
        if (assembled_ringbuf[i]) {
            // Here we make a copy of the shared_ptr, thus preserving the chunk
            ring.push_back(assembled_ringbuf[i]);
        }
    }
    pthread_mutex_unlock(&this->lock);
    return ring;
}

void assembled_chunk_ringbuf::put_unassembled_packet(const intensity_packet &packet, int64_t *event_counts)
{
    // We test these pointers instead of 'doneflag' so that we don't need to acquire the lock in every call.
    if (_unlikely(!active_chunk0 || !active_chunk1))
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");

    uint64_t packet_t0 = packet.fpga_count / packet.fpga_counts_per_sample;
    uint64_t packet_ichunk = packet_t0 / constants::nt_per_assembled_chunk;
    uint64_t active_ichunk = assembled_ringbuf_pos + assembled_ringbuf_size;

    if (packet_ichunk >= active_ichunk + 2) {
	//
	// If we receive a packet whose timestamps extend past the range of our current
	// assembly buffer, then we advance the buffer and send an assembled_chunk to the
	// "downstream" thread.
	//
	// A design decision here: for a packet which is far in the future, we advance the 
	// buffer by one assembled_chunk, rather than advancing all the way to the packet
	// timestamp.  This is to avoid a situation where a single rogue packet timestamped
	// in the far future effectively kills the L1 node.
	//
	this->_put_assembled_chunk(active_chunk0, event_counts);
	active_chunk0 = active_chunk1;
	active_chunk1 = this->_make_assembled_chunk(active_ichunk+2);
	active_ichunk++;
    }

    if (packet_ichunk == active_ichunk) {
	event_counts[intensity_network_stream::event_type::assembler_hit]++;
	active_chunk0->add_packet(packet);
    }
    else if (packet_ichunk == active_ichunk+1) {
	event_counts[intensity_network_stream::event_type::assembler_hit]++;
	active_chunk1->add_packet(packet);
    }
    else {
	event_counts[intensity_network_stream::event_type::assembler_miss]++;
	if (_unlikely(ini_params.throw_exception_on_assembler_miss))
	    throw runtime_error("ch_frb_io: assembler miss occurred, and this stream was constructed with the 'throw_exception_on_assembler_miss' flag");
    }
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

    if (assembled_ringbuf_size < constants::assembled_ringbuf_capacity) {
	// Add chunk to ring buffer
	int i = (assembled_ringbuf_pos + assembled_ringbuf_size) % constants::assembled_ringbuf_capacity;
	this->assembled_ringbuf[i] = chunk;
	this->assembled_ringbuf_size++;
	
	pthread_cond_broadcast(&this->cond_assembled_chunks_added);
	pthread_mutex_unlock(&this->lock);
	event_counts[intensity_network_stream::event_type::assembled_chunk_queued]++;
	return;
    }

    // If we get here, the ring buffer was full.
    pthread_mutex_unlock(&this->lock);
    event_counts[intensity_network_stream::event_type::assembled_chunk_dropped]++;

    if (ini_params.emit_warning_on_buffer_drop)
	cerr << "ch_frb_io: warning: processing thread is running too slow, dropping assembled_chunk\n";
    if (ini_params.throw_exception_on_buffer_drop)
	throw runtime_error("ch_frb_io: assembled_chunk was dropped and stream was constructed with 'throw_exception_on_buffer_drop' flag");
}


shared_ptr<assembled_chunk> assembled_chunk_ringbuf::get_assembled_chunk()
{
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (assembled_ringbuf_size > 0) {
	    int i = assembled_ringbuf_pos % constants::assembled_ringbuf_capacity;
	    shared_ptr<assembled_chunk> chunk = assembled_ringbuf[i];

	    this->assembled_ringbuf_pos++;
	    this->assembled_ringbuf_size--;
	    pthread_mutex_unlock(&this->lock);

	    if (!chunk)
		throw runtime_error("ch_frb_io: internal error: unexpected empty pointer in get_assembled_chunk()");

	    return chunk;
	}

	if (this->doneflag) {
	    // Ring buffer is empty and end_stream() has been called.
	    pthread_mutex_unlock(&this->lock);
	    return shared_ptr<assembled_chunk>();
	}
	
	// Wait for chunks to be added to the ring buffer.
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

    // Wake up processing thread, if it is waiting for data
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);

    this->doneflag = true;
    pthread_mutex_unlock(&this->lock);
}


std::shared_ptr<assembled_chunk> assembled_chunk_ringbuf::_make_assembled_chunk(uint64_t ichunk)
{
    if (ini_params.mandate_fast_kernels)
	return make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
    else if (ini_params.mandate_reference_kernels)
	return make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
    else
	return assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
}


}  // namespace ch_frb_io
