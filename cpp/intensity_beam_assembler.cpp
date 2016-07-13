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
    
    // wait for assembler thread to start
    pthread_mutex_lock(&ret->lock);
    while (!ret->assembler_thread_registered)
	pthread_cond_wait(&ret->cond_assembler_state_changed, &ret->lock);

    pthread_mutex_unlock(&ret->lock);    
    return ret;
}


intensity_beam_assembler::intensity_beam_assembler(int beam_id_) : beam_id(beam_id_)
{
    if ((beam_id < 0) || (beam_id >= 65536))
	throw runtime_error("intensity_beam_constructor: invalid beam_id");

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembler_state_changed, NULL);
    pthread_cond_init(&this->cond_unassembled_packets_added, NULL);
    pthread_cond_init(&this->cond_unassembled_packets_removed, NULL);

    for (int i = 0; i < unassembled_ringbuf_capacity; i++)
	unassembled_ringbuf[i].initialize(beam_id);
}


intensity_beam_assembler::~intensity_beam_assembler()
{
    for (int i = 0; i < unassembled_ringbuf_capacity; i++)
	unassembled_ringbuf[i].destroy();

    pthread_cond_destroy(&this->cond_assembler_state_changed);
    pthread_cond_destroy(&this->cond_unassembled_packets_added);
    pthread_cond_destroy(&this->cond_unassembled_packets_removed);
    pthread_mutex_destroy(&this->lock);
}


void intensity_beam_assembler::assembler_thread_register()
{
    pthread_mutex_lock(&this->lock);
    this->assembler_thread_registered = true;
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_mutex_unlock(&this->lock);
}


void intensity_beam_assembler::assembler_thread_unregister()
{
    pthread_mutex_lock(&this->lock);
    this->assembler_ended = true;
    this->assembler_thread_unregistered = true;
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_mutex_unlock(&this->lock);
}


void intensity_beam_assembler::end_assembler(bool join_thread)
{
    bool call_join_after_releasing_lock = false;

    pthread_mutex_lock(&this->lock);

    if (!assembler_thread_registered) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: intensity_beam_assembler::end_assembler() was called, but assembler_thread_registered flag wasn't set");
    }

    this->assembler_ended = true;

    // wake up all threads
    pthread_cond_broadcast(&this->cond_assembler_state_changed);
    pthread_cond_broadcast(&this->cond_unassembled_packets_added);
    pthread_cond_broadcast(&this->cond_unassembled_packets_removed);

    while (!assembler_thread_unregistered)
	pthread_cond_wait(&this->cond_assembler_state_changed, &this->lock);

    if (join_thread && !this->assembler_thread_joined) {
	this->assembler_thread_joined = true;
	call_join_after_releasing_lock = true;
    }

    pthread_mutex_unlock(&this->lock);

    if (call_join_after_releasing_lock)
	pthread_join(this->assembler_thread, NULL);
}


void intensity_beam_assembler::put_unassembled_packets(udp_packet_list &packet_list, bool blocking)
{
    if (packet_list.beam_id != this->beam_id)
	throw runtime_error("ch_frb_io: beam_id mismatch in intensity_beam_assembler::put_unassembled_packets()");

    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (assembler_ended) {
	    pthread_mutex_unlock(&this->lock);
	    throw runtime_error("ch_frb_io: intensity_beam_assembler::put_unassembled_packets() called, but assembler thread ended");
	}
	
	if (unassembled_ringbuf_size < unassembled_ringbuf_capacity)
	    break;

	if (!blocking) {
	    pthread_mutex_unlock(&this->lock);
	    cerr << "ch_frb_io: warning: assembler thread is running too slow, some packets will be dropped\n";
	    packet_list.clear();
	    return;
	}
	
	pthread_cond_wait(&this->cond_unassembled_packets_removed, &this->lock);
    }
	
    // success!
    int i = (unassembled_ringbuf_pos + unassembled_ringbuf_size) % unassembled_ringbuf_capacity;
    this->unassembled_ringbuf[i].swap(packet_list);
    this->unassembled_ringbuf_size++;

    pthread_cond_broadcast(&this->cond_unassembled_packets_added);
    pthread_mutex_unlock(&this->lock);
    packet_list.clear();
}


bool intensity_beam_assembler::get_unassembled_packets(udp_packet_list &packet_list)
{
    packet_list.clear();
    pthread_mutex_lock(&this->lock);

    while (unassembled_ringbuf_size == 0) {
	if (assembler_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_unassembled_packets_added, &this->lock);
    }

    int i = unassembled_ringbuf_pos % unassembled_ringbuf_capacity;
    this->unassembled_ringbuf[i].swap(packet_list);
    this->unassembled_ringbuf_pos++;
    this->unassembled_ringbuf_size--;

    pthread_cond_broadcast(&this->cond_unassembled_packets_removed);
    pthread_mutex_unlock(&this->lock);
    return true;
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

    assembler->assembler_thread_register();

    try {
	udp_packet_list packet_list;
	packet_list.initialize(assembler->beam_id);
	assembler_thread_main2(assembler.get(), packet_list);
	packet_list.destroy();
    } catch (...) {
	assembler->assembler_thread_unregister();
	throw;
    }

    assembler->assembler_thread_unregister();
    return NULL;
}


static void assembler_thread_main2(intensity_beam_assembler *assembler, udp_packet_list &packet_list)
{
    for (;;) {
	bool alive = assembler->get_unassembled_packets(packet_list);
	if (!alive)
	    return;

	cerr << "!";
    }
}


}  // namespace ch_frb_io
