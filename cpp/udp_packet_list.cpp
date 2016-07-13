#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


void udp_packet_list::initialize(int beam_id_)
{
    if ((beam_id_ < 0) || (beam_id_ >= 65536))
	throw runtime_error("ch_frb_io: udp_packet_list::initialize(): invalid beam_id");

    this->beam_id = beam_id_;
    this->npackets = 0;
    this->nbytes = 0;
    this->is_full = false;

    this->data_start = aligned_alloc<uint8_t> (max_bytes);
    this->data_end = data_start;

    this->packet_offsets = aligned_alloc<int> (max_packets+1);
    this->packet_offsets[0] = 0;
}


void udp_packet_list::destroy()
{
    free(data_start);
    free(packet_offsets);

    this->data_start = this->data_end = nullptr;
    this->packet_offsets = nullptr;
}


void udp_packet_list::swap(udp_packet_list &p)
{
    std::swap(this->beam_id, p.beam_id);
    std::swap(this->npackets, p.npackets);
    std::swap(this->nbytes, p.nbytes);
    std::swap(this->is_full, p.is_full);
    std::swap(this->data_start, p.data_start);
    std::swap(this->data_end, p.data_end);
    std::swap(this->packet_offsets, p.packet_offsets);
}


void udp_packet_list::add_packet(int packet_nbytes)
{
    // decided to pay a few cpu cycles for an extra check here...
    bool error = is_full || (packet_nbytes <= 0) || (packet_nbytes > max_packet_size);
    if (_unlikely(error))
	throw runtime_error("ch_frb_io: bad call to udp_packet_list::add_packet()");

    this->npackets++;
    this->nbytes += packet_nbytes;
    this->is_full = (npackets >= max_packets) || (nbytes + max_packet_size >= max_bytes);
    this->data_end = data_start + nbytes;
    this->packet_offsets[npackets] = nbytes;
}


void udp_packet_list::clear()
{
    this->npackets = 0;
    this->nbytes = 0;
    this->is_full = false;
    this->data_end = data_start;
    this->packet_offsets[0] = 0;
}


}  // ch_frb_io
