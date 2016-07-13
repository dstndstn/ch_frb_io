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
    while (!ret->assembler_started)
	pthread_cond_wait(&ret->cond_assembler_thread_started, &ret->lock);

    pthread_mutex_unlock(&ret->lock);    
    return ret;
}


// -------------------------------------------------------------------------------------------------


intensity_beam_assembler::intensity_beam_assembler(int beam_id_) : beam_id(beam_id_)
{
    if ((beam_id < 0) || (beam_id >= 65536))
	throw runtime_error("intensity_beam_constructor: invalid beam_id");

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembler_thread_started, NULL);
    pthread_cond_init(&this->cond_unassembled_packets_added, NULL);
    pthread_cond_init(&this->cond_unassembled_packets_removed, NULL);

    for (int i = 0; i < unassembled_ringbuf_capacity; i++)
	unassembled_ringbuf[i].initialize(beam_id);
}


intensity_beam_assembler::~intensity_beam_assembler()
{
    for (int i = 0; i < unassembled_ringbuf_capacity; i++)
	unassembled_ringbuf[i].destroy();

    pthread_cond_destroy(&this->cond_assembler_thread_started);
    pthread_cond_destroy(&this->cond_unassembled_packets_added);
    pthread_cond_destroy(&this->cond_unassembled_packets_removed);
    pthread_mutex_destroy(&this->lock);
}


void intensity_beam_assembler::end_assembler(bool join_thread)
{
    bool call_join_after_releasing_lock = false;

    pthread_mutex_lock(&this->lock);
    this->assembler_ended = true;

    if (join_thread && !this->assembler_joined) {
	this->assembler_joined = true;
	call_join_after_releasing_lock = true;
    }

    // wake up any waiting threads
    pthread_cond_broadcast(&this->cond_unassembled_packets_added);
    pthread_cond_broadcast(&this->cond_unassembled_packets_removed);
    pthread_mutex_unlock(&this->lock);

    if (call_join_after_releasing_lock)
	pthread_join(this->assembler_thread, NULL);
}


void intensity_beam_assembler::put_unassembled_packets(udp_packet_list &packet_list_to_swap, bool blocking)
{
    if (packet_list_to_swap.beam_id != this->beam_id)
	throw runtime_error("ch_frb_io: beam_id mismatch in intensity_beam_assembler::put_unassembled_packets()");

    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (assembler_ended) {
	    pthread_mutex_unlock(&this->lock);
	    throw runtime_error("intensity_beam_assembler::put_unassembled_packets() called, but assembler thread ended");
	}
	
	if (unassembled_ringbuf_size < unassembled_ringbuf_capacity)
	    break;

	if (!blocking) {
	    pthread_mutex_unlock(&this->lock);
	    cerr << "ch_frb_io: warning: assembler thread is running too slow, some packets will be dropped\n";
	    packet_list_to_swap.clear();
	    return;
	}
	
	pthread_cond_wait(&this->cond_unassembled_packets_removed, &this->lock);
    }
	
    // success!
    int i = (unassembled_ringbuf_pos + unassembled_ringbuf_size) % unassembled_ringbuf_capacity;
    this->unassembled_ringbuf[i].swap(packet_list_to_swap);
    this->unassembled_ringbuf_size++;

    pthread_cond_broadcast(&this->cond_unassembled_packets_added);
    pthread_mutex_unlock(&this->lock);
    packet_list_to_swap.clear();
}


bool intensity_beam_assembler::get_unassembled_packets(udp_packet_list &packet_list_to_swap)
{
    packet_list_to_swap.clear();
    pthread_mutex_lock(&this->lock);

    while (unassembled_ringbuf_size == 0) {
	if (assembler_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_unassembled_packets_added, &this->lock);
    }

    // return packets to caller
    int i = unassembled_ringbuf_pos % unassembled_ringbuf_capacity;
    this->unassembled_ringbuf[i].swap(packet_list_to_swap);
    this->unassembled_ringbuf_pos++;
    this->unassembled_ringbuf_size--;

    pthread_mutex_unlock(&this->lock);
    return true;
}


void intensity_beam_assembler::assembler_thread_register()
{
    pthread_mutex_lock(&this->lock);
    this->assembler_started = true;
    pthread_cond_broadcast(&this->cond_assembler_thread_started);
    pthread_mutex_unlock(&this->lock);
}


// -------------------------------------------------------------------------------------------------


static void *assembler_thread_main(void *opaque_arg)
{
    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<intensity_beam_assembler> *arg = (shared_ptr<intensity_beam_assembler> *) opaque_arg;
    shared_ptr<intensity_beam_assembler> assembler = *arg;   // 'exchanger' is safe to use below
    delete arg;

    assembler->assembler_thread_register();

    udp_packet_list packet_list;
    packet_list.initialize(assembler->beam_id);

    for (;;) {
	bool alive = assembler->get_unassembled_packets(packet_list);
	if (!alive) {
	    packet_list.destroy();
	    return NULL;
	}

	cerr << "!!! assembler thread got packet list: [";
	for (int i = 0; i < packet_list.npackets; i++)
	    cerr << " " << (packet_list.packet_offsets[i+1] - packet_list.packet_offsets[i]);
	cerr << " ]\n";
    }
}


}  // namespace ch_frb_io
