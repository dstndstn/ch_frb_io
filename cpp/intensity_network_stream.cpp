// General FIXME: currently, we silently drop bad packets

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// FIXME also in intensity_network_ostream.cpp
inline int packet_size(int nbeam, int nfreq, int nupfreq, int ntsamp)
{
    return 24 + 2*nbeam + 2*nfreq + 8*nbeam*nfreq + (nbeam * nfreq * nupfreq * ntsamp);
}


// -------------------------------------------------------------------------------------------------
//
// L0L1_packet_list


void L0L1_packet_list::initialize(int beam_id_)
{
    if ((beam_id_ < 0) || (beam_id_ >= 65536))
	throw runtime_error("ch_frb_io: L0L1_packet_list::initialize(): invalid beam_id");

    this->beam_id = beam_id_;
    this->npackets = 0;
    this->nbytes = 0;
    this->data_buf = aligned_alloc<uint8_t> (max_bytes);
    this->packet_offsets = aligned_alloc<int> (max_packets+1);
    this->packet_offsets[0] = 0;
}


void L0L1_packet_list::destroy()
{
    free(data_buf);
    free(packet_offsets);

    this->npackets = this->nbytes = 0;
    this->packet_offsets = nullptr;
    this->data_buf = nullptr;
}


bool L0L1_packet_list::add_packet(int packet_nbytes)
{
    // decided to pay a few cpu cycles for an extra check here...
    bool internal_error = ((packet_nbytes <= 0) || 
			   (packet_nbytes > max_packet_size) ||
			   (this->npackets >= max_packets) || 
			   (this->nbytes >= (max_bytes - max_packet_size)));

    if (_unlikely(internal_error))
	throw runtime_error("ch_frb_io: internal error in L0L1_packet_list::add_packet()");

    this->npackets++;
    this->nbytes += packet_nbytes;
    this->packet_offsets[npackets] = nbytes;

    bool is_full = (npackets >= max_packets) || (nbytes >= (max_bytes - max_packet_size));
    return is_full;
}


void L0L1_packet_list::clear()
{
    this->npackets = 0;
    this->nbytes = 0;
    this->packet_offsets[0] = 0;
}


// -------------------------------------------------------------------------------------------------
//
// intensity_beam_assembler


intensity_beam_assembler::intensity_beam_assembler(int beam_id_) :
    beam_id(beam_id_)
{ 
    if ((beam_id < 0) || (beam_id >= 65536))
	throw runtime_error("ch_frb_io: intensity_beam_assembler constructor: invalid beam_id");
}


void intensity_beam_assembler::send_packet_list(L0L1_packet_list &packet_list)
{
    if (packet_list.beam_id != this->beam_id)
	throw runtime_error("ch_frb_io: intensity_beam_assembler::send_packet_list(): beam_id mismatch");

    // placeholder
    cerr << "!!! assembler received packet list, npackets=" << packet_list.npackets << " !!!\n";
    packet_list.npackets = 0;
    packet_list.nbytes = 0;
}


// -------------------------------------------------------------------------------------------------
//
// intensity_packet_stream


intensity_packet_stream::~intensity_packet_stream()
{
    for (int i = 0; i < nassemblers; i++)
	this->packet_lists[i].destroy();

    this->beam_assemblers.clear();
    this->nassemblers = 0;
}


void intensity_packet_stream::add_beam(const shared_ptr<intensity_beam_assembler> &b)
{
    if (network_thread_valid)
	throw runtime_error("ch_frb_io: add_beam() called on packet stream with network thread already running");
    if (nassemblers >= max_beams)
	throw runtime_error("ch_frb_io: intensity_packet_stream::add_beam(): max number of beams has already been reached");

    int beam_id = b->beam_id;
    for (int i = 0; i < nassemblers; i++)
	if (packet_lists[i].beam_id == beam_id)
	    throw runtime_error("ch_frb_io: intensity_packet_streram::add_beam(): duplicate beam_id");
    
    this->packet_lists[nassemblers].initialize(beam_id);
    this->beam_assemblers.push_back(b);
    this->nassemblers++;
}


bool intensity_packet_stream::process_packet(int packet_nbytes, const uint8_t *in)
{
    if (_unlikely(nassemblers == 0))
	throw runtime_error("ch_frb_io: intensity_packet_stream::process_packet() called on stream with no assemblers");

    if (_unlikely(packet_nbytes < 24))
	return true;  // silently drop bad packets

    uint32_t protocol_version = *((uint32_t *) in);
    if (_unlikely(protocol_version != 1))
	return true;  // silently drop bad packets

    int data_nbytes = *((int16_t *) (in+4));
    int packet_nbeam = *((uint16_t *) (in+16));
    int packet_nfreq = *((uint16_t *) (in+18));
    int packet_nupfreq = *((uint16_t *) (in+20));
    int packet_ntsamp = *((uint16_t *) (in+22));

    if (packet_nbytes == 24)
	return false;   // this is the only path which returns false (indicating end of stream)

    // The following way of writing the comparisons (using uint64_t) guarantees no overflow
    uint64_t n4 = uint64_t(packet_nbeam) * uint64_t(packet_nfreq) * uint64_t(packet_nupfreq) * uint64_t(packet_ntsamp);
    if (_unlikely(!n4 || (n4 > 10000)))
	return true;    // silently drop bad packets
    if (_unlikely(data_nbytes != int(n4)))
	return true;    // silently drop bad packets
    if (_unlikely(packet_nbytes != packet_size(packet_nbeam, packet_nfreq, packet_nupfreq, packet_ntsamp)))
	return true;    // silently drop bad packets

    const uint16_t *packet_beam_ids = (const uint16_t *) (in + 24);
    const uint8_t *packet_freq_ids = in + 24 + 2*packet_nbeam;
    const uint8_t *packet_scales = in + 24 + 2*packet_nbeam + 2*packet_nfreq;
    const uint8_t *packet_offsets = in + 24 + 2*packet_nbeam + 2*packet_nfreq + 4*packet_nbeam*packet_nfreq;
    const uint8_t *packet_data = in + 24 + 2*packet_nbeam + 2*packet_nfreq + 8*packet_nbeam*packet_nfreq;
    
    // Loop over beams in the packet, matching to beam_assembler objects.

    for (int ibeam = 0; ibeam < packet_nbeam; ibeam++) {
	for (int jbeam = 0; jbeam < ibeam; jbeam++)
	    if (_unlikely(packet_beam_ids[ibeam] == packet_beam_ids[jbeam]))
		return true;   // duplicate beam_ids in packet, silently drop
	
	// Find assembler matching beam ID

	int beam_id = packet_beam_ids[ibeam];
	int assembler_ix = 0;

	for (;;) {
	    if (assembler_ix >= nassemblers)
		return true;   // bad beam_id in packet, silently drop
	    if (packet_lists[assembler_ix].beam_id == beam_id)
		break;          // found assembler
	    assembler_ix++;
	}
	
	// Make subpacket
	
	uint8_t *subpacket = packet_lists[assembler_ix].data_buf + packet_lists[assembler_ix].nbytes;
	uint8_t *subpacket_freq_ids = subpacket + 26;
	uint8_t *subpacket_scales = subpacket + 26 + 2*packet_nfreq;
	uint8_t *subpacket_offsets = subpacket + 26 + 6*packet_nfreq;
	uint8_t *subpacket_data = subpacket + 26 + 10*packet_nfreq;

	int subpacket_data_size = packet_nfreq * packet_nupfreq * packet_ntsamp;
	int subpacket_total_size = subpacket_data_size + 26 + 10*packet_nfreq;

	memcpy(subpacket, in, 24);                                // copy header and overwrite a few fields...
	*((int16_t *) (subpacket+4)) = subpacket_data_size;       // overwrite 'data_nbytes' field
	*((uint16_t *) (subpacket+16)) = uint16_t(1);             // overwrite 'nbeams' field
	*((uint16_t *) (subpacket+24)) = packet_beam_ids[ibeam];  // beam_id field

	memcpy(subpacket_freq_ids, packet_freq_ids, 2*packet_nfreq);
	memcpy(subpacket_scales, packet_scales + 4*ibeam*packet_nfreq, 4*packet_nfreq);
	memcpy(subpacket_offsets, packet_offsets + 4*ibeam*packet_nfreq, 4*packet_nfreq);
	memcpy(subpacket_data, packet_data + ibeam*subpacket_data_size, subpacket_data_size);

	// Record the subpacket in the packet_list
	bool packet_list_is_full = packet_lists[assembler_ix].add_packet(subpacket_total_size);
	
	// If packet_list is full, hand it off to assembler thread (nonblocking).
	// When send_packet_list() returns, the packet list will be empty and ready to receive new packets.
	if (packet_list_is_full)
	    beam_assemblers[assembler_ix]->send_packet_list(packet_lists[assembler_ix]);
    }

    return true;
}


void intensity_packet_stream::finalize()
{
    for (int i = 0; i < nassemblers; i++) {
	if (packet_lists[i].npackets > 0)
	    beam_assemblers[i]->send_packet_list(packet_lists[i]);
	packet_lists[i].destroy();
    }

    this->beam_assemblers.clear();
    this->nassemblers = 0;
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


struct network_thread_context {
    uint16_t udp_port;
    shared_ptr<intensity_packet_stream> stream;

    network_thread_context(uint16_t udp_port_, const shared_ptr<intensity_packet_stream> &stream_) :
	udp_port(udp_port_), stream(stream_)
    { }
};


static void *network_thread_main(void *opaque_arg)
{
    // FIXME is 2MB socket_bufsize a good choice?  I don't see anything wrong with using a larger value, but 2MB is the max allowed on my osx laptop!
    static constexpr ssize_t max_packet_size = 16384;
    static const int socket_bufsize = 2 << 21;   // 2 MB

    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_thread_main()");

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<network_thread_context> *arg = (shared_ptr<network_thread_context> *) opaque_arg;
    shared_ptr<network_thread_context> context = *arg;   // 'context' is safe to use below
    delete arg;

    if (!context)
	throw runtime_error("ch_frb_io: internal error: no network_thread_context passed to network_thread_main()");
    if (!context->stream)
	throw runtime_error("ch_frb_io: internal error: no intensity_packet_stream passed to network_thread_main()");    

    // Bare pointer (note that reference is held through context->stream shared_ptr)
    intensity_packet_stream *stream = context->stream.get();

    int sock_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (sock_fd < 0)
	throw runtime_error(string("ch_frb_io: socket() failed: ") + strerror(errno));

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));
	
    server_address.sin_family = AF_INET;
    inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
    server_address.sin_port = htons(context->udp_port);
    
    int err = ::bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: bind() failed: ") + strerror(errno));

    err = setsockopt(sock_fd, SOL_SOCKET, SO_RCVBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt() failed: ") + strerror(errno));

    vector<uint8_t> packet_buf(max_packet_size);

    cerr << ("ch_frb_io: network thread listening for packets on port " + to_string(context->udp_port) + "\n");
    
    for (;;) {
	ssize_t bytes_read = read(sock_fd, (char *) &packet_buf[0], max_packet_size);

	if (bytes_read >= max_packet_size)
	    continue;  // silently drop bad packet

	bool end_of_stream = !stream->process_packet(bytes_read, &packet_buf[0]);
	if (end_of_stream)
	    break;
    }

    cerr << "ch_frb_io: end of stream, network thread exiting\n";

    stream->finalize();
    return NULL;
}


void spawn_network_thread(int udp_port, const shared_ptr<intensity_packet_stream> &stream)
{
    if (stream->network_thread_valid)
	throw runtime_error("ch_frb_io: spawn_network_thread() called on packet stream with network thread already running");
    if ((udp_port <= 0) || (udp_port >= 65536))
	throw runtime_error("ch_frb_io: spawn_network_thread(): bad udp port " + to_string(udp_port));

    shared_ptr<network_thread_context> context = make_shared<network_thread_context> (udp_port, stream);

    // To hand off the context to the network thread, we use a bare pointer to a shared pointer
    shared_ptr<network_thread_context> *p = new shared_ptr<network_thread_context> (context);

    int err = pthread_create(&stream->network_thread, NULL, network_thread_main, p);
    
    if (err) {
	delete p;   // thread didn't start, caller is responsible for deleting p
	throw runtime_error("ch_frb_io: spawn_network_thread(): pthread_create failed");
    }

    // network thread started and is responsible for deleting p
    stream->network_thread_valid = true;
}


void wait_for_end_of_stream(const shared_ptr<intensity_packet_stream> &stream)
{
    if (!stream->network_thread_valid)
	throw runtime_error("ch_frb_io: wait_for_end_of_stream() called on packet stream with no network thread");

    pthread_join(stream->network_thread, NULL);
    stream->network_thread_valid = false;
}


}  // namespace ch_frb_io
