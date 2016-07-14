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
static void  assembler_thread_main2(intensity_beam_assembler *assembler, udp_packet_list &packet_list);


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


intensity_beam_assembler::intensity_beam_assembler(int beam_id_) : beam_id(beam_id_)
{
    if ((beam_id < 0) || (beam_id >= 65536))
	throw runtime_error("intensity_beam_constructor: invalid beam_id");

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembler_state_changed, NULL);
    pthread_cond_init(&this->cond_unassembled_packets_added, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);

    for (int i = 0; i < unassembled_ringbuf_capacity; i++)
	unassembled_ringbuf[i].initialize(beam_id);
}


intensity_beam_assembler::~intensity_beam_assembler()
{
    for (int i = 0; i < unassembled_ringbuf_capacity; i++)
	unassembled_ringbuf[i].destroy();

    pthread_cond_destroy(&this->cond_assembler_state_changed);
    pthread_cond_destroy(&this->cond_unassembled_packets_added);
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
    pthread_cond_broadcast(&this->cond_unassembled_packets_added);
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);

    if (join_thread && !this->assembler_thread_joined) {
	this->assembler_thread_joined = true;
	call_join_after_releasing_lock = true;
    }

    pthread_mutex_unlock(&this->lock);

    if (call_join_after_releasing_lock)
	pthread_join(this->assembler_thread, NULL);
}


// XXX don't forget to write comment here explaining semantics of the 'packet_list' arg
bool intensity_beam_assembler::put_unassembled_packets(udp_packet_list &packet_list)
{
    if (_unlikely(packet_list.beam_id != this->beam_id))
	throw runtime_error("ch_frb_io: internal error: beam_id mismatch in intensity_beam_assembler::put_unassembled_packets()");

    pthread_mutex_lock(&this->lock);

    if (_unlikely(!stream_started)) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: intensity_beam_assembler::put_unassembled_packets() was called, but start_stream() was never called");
    }

    if (stream_ended) {
	pthread_mutex_unlock(&this->lock);
	return false;
    }
	
    if (unassembled_ringbuf_size >= unassembled_ringbuf_capacity) {
	pthread_mutex_unlock(&this->lock);
	cerr << "ch_frb_io: warning: assembler thread is running too slow, some packets will be dropped\n";
	packet_list.clear();
	return true;
    }
	
    int i = (unassembled_ringbuf_pos + unassembled_ringbuf_size) % unassembled_ringbuf_capacity;
    this->unassembled_ringbuf[i].swap(packet_list);
    this->unassembled_ringbuf_size++;

    pthread_cond_broadcast(&this->cond_unassembled_packets_added);
    pthread_mutex_unlock(&this->lock);
    packet_list.clear();
    return true;
}


bool intensity_beam_assembler::get_unassembled_packets(udp_packet_list &packet_list)
{
    packet_list.clear();
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (unassembled_ringbuf_size > 0) {
	    int i = unassembled_ringbuf_pos % unassembled_ringbuf_capacity;
	    this->unassembled_ringbuf[i].swap(packet_list);
	    this->unassembled_ringbuf_pos++;
	    this->unassembled_ringbuf_size--;
	    
	    pthread_mutex_unlock(&this->lock);
	    return true;
	}

	if (stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_unassembled_packets_added, &this->lock);
    }
}


void intensity_beam_assembler::put_assembled_chunk(const shared_ptr<assembled_chunk> &chunk)
{
    pthread_mutex_lock(&this->lock);

    if (assembled_ringbuf_size >= assembled_ringbuf_capacity) {
	pthread_mutex_unlock(&this->lock);
	cerr << "ch_frb_io: warning: post-assembler thread is running too slow, some packets will be dropped\n";
	return;
    }

    int i = (assembled_ringbuf_pos + assembled_ringbuf_size) % assembled_ringbuf_capacity;
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
	    chunk = assembled_ringbuf[assembled_ringbuf_pos % assembled_ringbuf_capacity];
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

    assembler->assembler_thread_startup();

    try {
	udp_packet_list packet_list;
	packet_list.initialize(assembler->beam_id);
	assembler_thread_main2(assembler.get(), packet_list);
	packet_list.destroy();
    } catch (...) {
	assembler->end_stream(false);  // the "false" means "nonblocking" (otherwise we'd deadlock!)
	throw;
    }

    assembler->end_stream(false);  // "false" means same as above
    return NULL;
}


static void assembler_thread_main2(intensity_beam_assembler *assembler, udp_packet_list &unassembled_packet_list)
{
    int fpga_counts_per_sample, nupfreq;
    assembler->wait_for_stream_params(fpga_counts_per_sample, nupfreq);
    
    // Outer loop over unassembled packets

    for (;;) {
	bool alive = assembler->get_unassembled_packets(unassembled_packet_list);
	if (!alive)
	    return;

	cerr << "!";
    }
}


}  // namespace ch_frb_io
