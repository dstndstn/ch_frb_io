#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


static_assert(constants::max_input_udp_packet_size >= constants::max_output_udp_packet_size,
	      "udp_packet_list.cpp assumes constants::max_input_udp_packet_size >= constants::max_output_udp_packet_size");


// -------------------------------------------------------------------------------------------------
//
// udp_packet_list


udp_packet_list::udp_packet_list(int max_npackets_, int max_nbytes_) :
    max_npackets(max_npackets_), max_nbytes(max_nbytes_)
{
    if (max_npackets <= 0)
	throw runtime_error("udp_packet_list constructor: expected max_npackets > 0");
    if (max_nbytes <= 0)
	throw runtime_error("udp_packet_list constructor: expected max_nbytes > 0");

    // Note: this->curr_npackets and this->curr_nbytes are initialized to zero automatically.
    this->buf = unique_ptr<uint8_t[]> (new uint8_t[max_nbytes + constants::max_input_udp_packet_size]);
    this->off_buf = unique_ptr<int[]> (new int[max_npackets + 1]);
    this->is_full = false;
    this->data_start = buf.get();
    this->data_end = data_start;
    this->packet_offsets = off_buf.get();
    this->packet_offsets[0] = 0;
}


void udp_packet_list::add_packet(int packet_nbytes)
{
    if ((packet_nbytes <= 0) || (packet_nbytes > constants::max_input_udp_packet_size))
	throw runtime_error("udp_packet_list::add_packet(): bad value of 'packet_nbytes'");
    if (is_full)
	throw runtime_error("udp_packet_list::add_packet() called on full packet_list");

    this->curr_npackets++;
    this->curr_nbytes += packet_nbytes;
    this->is_full = (curr_npackets >= max_npackets) || (curr_nbytes >= max_nbytes);
    this->data_end = data_start + curr_nbytes;
    this->packet_offsets[curr_npackets] = curr_nbytes;
}


void udp_packet_list::reset()
{
    this->curr_npackets = 0;
    this->curr_nbytes = 0;
    this->is_full = false;
    this->data_end = data_start;
    this->packet_offsets[0] = 0;
}


// -------------------------------------------------------------------------------------------------
//
// udp_packet_ringbuf


udp_packet_ringbuf::udp_packet_ringbuf(int ringbuf_capacity_, int max_npackets_per_list_, int max_nbytes_per_list_, const std::string &dropmsg_, bool drops_allowed_)
    : drops_allowed(drops_allowed_),
      ringbuf_capacity(ringbuf_capacity_), 
      max_npackets_per_list(max_npackets_per_list_),
      max_nbytes_per_list(max_nbytes_per_list_), 
      dropmsg(dropmsg_)
{
    if (ringbuf_capacity <= 0)
	throw runtime_error("udp_packet_ringbuf constructor: expected ringbuf_capacity > 0");

    this->ringbuf.resize(ringbuf_capacity);

    for (int i = 0; i < ringbuf_capacity; i++) {
	udp_packet_list l(max_npackets_per_list, max_nbytes_per_list);
	std::swap(this->ringbuf[i], l);
    }

    xpthread_mutex_init(&this->lock);
    xpthread_cond_init(&this->cond_packets_added);
    xpthread_cond_init(&this->cond_packets_removed);
}


udp_packet_ringbuf::~udp_packet_ringbuf()
{
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&cond_packets_added);
    pthread_cond_destroy(&cond_packets_removed);
}


bool udp_packet_ringbuf::producer_put_packet_list(udp_packet_list &packet_list, bool is_blocking)
{    
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	if (ringbuf_size < ringbuf_capacity) {
	    int i = (ringbuf_pos + ringbuf_size) % ringbuf_capacity;
	    std::swap(this->ringbuf[i], packet_list);
	    this->ringbuf_size++;
	
	    pthread_cond_broadcast(&this->cond_packets_added);
	    pthread_mutex_unlock(&this->lock);
	    packet_list.reset();
	    return true;
	}

	if (!is_blocking) {
	    pthread_mutex_unlock(&this->lock);
	    packet_list.reset();
	    
	    if (dropmsg.size() > 0)
		cerr << dropmsg << endl;

	    if (!drops_allowed) {
		cerr << "ring buffer overfilled, and 'drops_allowed' flag was set to false, this will be treated as an error!\n";
		exit(1);
	    }
	    
	    return true;
	}

	pthread_cond_wait(&this->cond_packets_removed, &this->lock);
    }
}


bool udp_packet_ringbuf::consumer_get_packet_list(udp_packet_list &packet_list)
{
    packet_list.reset();
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (ringbuf_size > 0) {
	    int i = ringbuf_pos % ringbuf_capacity;
	    std::swap(this->ringbuf[i], packet_list);
	    this->ringbuf_pos++;
	    this->ringbuf_size--;

	    pthread_cond_broadcast(&this->cond_packets_removed);
	    pthread_mutex_unlock(&this->lock);
	    return true;
	}

	if (stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_packets_added, &this->lock);
    }
}


void udp_packet_ringbuf::end_stream()
{
    pthread_mutex_lock(&this->lock);
    this->stream_ended = true;
    pthread_cond_broadcast(&this->cond_packets_added);
    pthread_cond_broadcast(&this->cond_packets_removed);
    pthread_mutex_unlock(&this->lock);
}


udp_packet_list udp_packet_ringbuf::allocate_packet_list() const
{
    return udp_packet_list(this->max_npackets_per_list, this->max_nbytes_per_list);
}


}  // namespace ch_frb_io
