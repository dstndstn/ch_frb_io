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


// Defined later in this file
static void *network_thread_main(void *opaque_arg);
static ssize_t network_thread_main2(intensity_network_stream *stream, int sock_fd);


// -------------------------------------------------------------------------------------------------
//
// intensity_network_stream


// Static member function (de facto constructor)
shared_ptr<intensity_network_stream> intensity_network_stream::make(const vector<shared_ptr<intensity_beam_assembler> > &assemblers, int udp_port)
{
    intensity_network_stream *ret_bare = new intensity_network_stream(assemblers, udp_port);
    shared_ptr<intensity_network_stream> ret(ret_bare);
    
    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<intensity_network_stream> *arg = new shared_ptr<intensity_network_stream> (ret);
    int err = pthread_create(&ret->network_thread, NULL, network_thread_main, arg);

    if (err) {
        delete arg;   // note: if pthread_create() succeeds, then assembler thread will delete this pointer
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_stream constructor: ") + strerror(errno));
    }
    
    // wait for network thread to start
    pthread_mutex_lock(&ret->lock);
    while (!ret->network_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->lock);

    pthread_mutex_unlock(&ret->lock);    
    return ret;
}


intensity_network_stream::intensity_network_stream(const vector<shared_ptr<intensity_beam_assembler> > &assemblers_, int udp_port_) :
    assemblers(assemblers_), udp_port(udp_port_)
{
    // Argument checking

    int nassemblers = assemblers.size();

    if (nassemblers == 0)
	throw runtime_error("ch_frb_io: empty assembler list passed to intensity_network_stream constructor");

    for (int i = 0; i < nassemblers; i++) {
	if (!assemblers[i])
	    throw runtime_error("ch_frb_io: empty assembler pointer passed to intensity_network_stream constructor");

	for (int j = 0; j < i; j++)
	    if (assemblers[i]->beam_id == assemblers[j]->beam_id)
		throw runtime_error("ch_frb_io: assembler list passed to intensity_network_stream constructor contains duplicate beam_ids");
    }

    if ((udp_port <= 0) || (udp_port >= 65536))
	throw runtime_error("ch_frb_io: intensity_network_stream constructor: bad udp port " + to_string(udp_port));

    // Initialize pthread members

    pthread_mutex_init(&lock, NULL);
    pthread_cond_init(&cond_state_changed, NULL);
}


intensity_network_stream::~intensity_network_stream()
{
    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&lock);
}


void intensity_network_stream::network_thread_startup()
{
    pthread_mutex_lock(&this->lock);
    this->network_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->lock);
}


void intensity_network_stream::wait_for_network_thread_startup()
{
    pthread_mutex_lock(&this->lock);

    while (!network_thread_started)
	pthread_cond_wait(&this->cond_state_changed, &this->lock);

    pthread_mutex_unlock(&this->lock);
}


void intensity_network_stream::start_stream()
{
    pthread_mutex_lock(&this->lock);

    if (stream_started) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::start_stream() called on running, completed, or cancelled stream");
    }

    this->stream_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->lock);
}


bool intensity_network_stream::wait_for_start_stream()
{
    pthread_mutex_lock(&this->lock);

    while (!stream_started)
	pthread_cond_wait(&this->cond_state_changed, &this->lock);

    bool retval = !stream_ended;
    pthread_mutex_unlock(&this->lock);
    return retval;
}


void intensity_network_stream::end_stream()
{
    pthread_mutex_lock(&this->lock);

    // Set flags as if stream had run to completion.  This is convenient e.g. for waking up threads
    // which are blocked in wait_for_packets().
    this->stream_started = true;
    this->stream_ended = true;
    
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->lock);

    for (unsigned int i = 0; i < assemblers.size(); i++)
	assemblers[i]->end_stream();
}


void intensity_network_stream::join_all_threads()
{
    bool call_join_after_releasing_lock = false;

    pthread_mutex_lock(&lock);
    
    if (!stream_started) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::wait_for_end_of_stream() was called with no prior call to start_stream()");
    }
    
    while (!stream_ended)
	pthread_cond_wait(&this->cond_state_changed, &this->lock);

    if (!network_thread_joined) {
	call_join_after_releasing_lock = true;
	this->network_thread_joined = true;
    }

    pthread_mutex_unlock(&lock);

    if (call_join_after_releasing_lock)
	pthread_join(network_thread, NULL);

    for (unsigned int i = 0; i < assemblers.size(); i++)
	assemblers[i]->join_assembler_thread();
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


static void *network_thread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_thread_main()");

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<intensity_network_stream> *arg = (shared_ptr<intensity_network_stream> *) opaque_arg;
    shared_ptr<intensity_network_stream> stream = *arg;
    delete arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to network_thread_main()");

    stream->network_thread_startup();

    int sock_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (sock_fd < 0)
	throw runtime_error(string("ch_frb_io: socket() failed: ") + strerror(errno));

    ssize_t npackets_received = 0;

    // We use a try..catch to ensure that the socket always gets closed, and end_stream()
    // always gets called, even if an exception is thrown.

    try {
	npackets_received = network_thread_main2(stream.get(), sock_fd);
    } catch (...) {
	close(sock_fd);
	stream->end_stream();
	throw;
    }

    close(sock_fd);
    stream->end_stream();

    return NULL;
}


// Returns number of packets received
static ssize_t network_thread_main2(intensity_network_stream *stream, int sock_fd)
{
    int nassemblers = stream->assemblers.size();

    // We unpack the assemblers into an array of bare pointers, and unpack the beam_ids
    // into an array of ints, for speed, although it probably doesn't matter!

    vector<intensity_beam_assembler *> assemblers(nassemblers, NULL);
    vector<int> assembler_beam_ids(nassemblers);

    for (int i = 0; i < nassemblers; i++) {
	assemblers[i] = stream->assemblers[i].get();
	assembler_beam_ids[i] = stream->assemblers[i]->beam_id;
    }

    // Socket setup: bind(), bufsize, timeout

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));
	
    server_address.sin_family = AF_INET;
    inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
    server_address.sin_port = htons(stream->udp_port);
    
    int err = ::bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: bind() failed: ") + strerror(errno));

    int socket_bufsize = constants::recv_socket_bufsize;
    err = setsockopt(sock_fd, SOL_SOCKET, SO_RCVBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVBUF) failed: ") + strerror(errno));

    struct timeval tv_timeout = { 0, constants::recv_socket_timeout_usec };
    err = setsockopt(sock_fd, SOL_SOCKET, SO_RCVTIMEO, &tv_timeout, sizeof(tv_timeout));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVTIME0) failed: ") + strerror(errno));

    cerr << ("ch_frb_io: network thread: listening for packets on port " + to_string(stream->udp_port) + "\n");

    // Wait for stream to start

    bool cancelled = !stream->wait_for_start_stream();
    
    if (cancelled)
	return 0;

    vector<udp_packet_list> assembler_packet_lists(nassemblers);
    for (int i = 0; i < nassemblers; i++)
	assembler_packet_lists[i] = assemblers[i]->allocate_unassembled_packet_list();

    struct timeval tv_ini = xgettimeofday();
    vector<int64_t> assembler_timestamps(nassemblers, 0);   // microseconds relative to tv_ini

    vector<uint8_t> packet_buf(constants::max_input_udp_packet_size + 1);
    uint8_t *packet = &packet_buf[0];
    ssize_t npackets_received = 0;
    
    //
    // Main packet loop!
    // 
    // FIXME there is a lot of silent "swallowing" of error conditions here, needs revisiting to
    // think about what really makes sense.
    //

    while (stream->is_alive()) {
	// Check whether any assembled_packet_lists have timed out.
	int64_t curr_timestamp = usec_between(tv_ini, xgettimeofday());
	int64_t threshold_timestamp = curr_timestamp - constants::unassembled_ringbuf_timeout_usec;

	for (int i = 0; i < nassemblers; i++) {
	    if (assembler_timestamps[i] > threshold_timestamp)
		continue;
	    if (assembler_packet_lists[i].curr_npackets == 0)
		continue;

	    bool assembler_alive = assemblers[i]->put_unassembled_packets(assembler_packet_lists[i]);
		
	    if (!assembler_alive) {
		cerr << "ch_frb_io:  assembler thread died unexpectedly!  network thread will die too...\n";
		return npackets_received;
	    }
	}

	// Read new packet from socket (note that read() can time out)
	int packet_nbytes = read(sock_fd, (char *) packet, constants::max_input_udp_packet_size + 1);

	if (packet_nbytes < 0) {
	    if ((errno == EAGAIN) || (errno == ETIMEDOUT))
		continue;  // timed out
	    throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
	}

	if (_unlikely(packet_nbytes > constants::max_input_udp_packet_size))
	    continue;  // silently drop bad packet

	uint32_t protocol_version = *((uint32_t *) packet);
	if (_unlikely(protocol_version != 1))
	    continue;  // silently drop bad packet

	int data_nbytes = *((int16_t *) (packet+4));
	int packet_fpga_counts_per_sample = *((uint16_t *) (packet+6));
	int packet_nbeam = *((uint16_t *) (packet+16));
	int packet_nfreq = *((uint16_t *) (packet+18));
	int packet_nupfreq = *((uint16_t *) (packet+20));
	int packet_ntsamp = *((uint16_t *) (packet+22));
	
	// If we receive a special "short" packet (length 24), it indicates end-of-stream.
	// FIXME is this a temporary kludge or something which should be documented in the packet protocol?
	if (_unlikely(packet_nbytes == 24)) {
	    cerr << "ch_frb_io: network thread received end-of-stream packets\n";
	    break;
	}
	
	// The following way of writing the comparisons (using uint64_t) guarantees no overflow
	uint64_t n4 = uint64_t(packet_nbeam) * uint64_t(packet_nfreq) * uint64_t(packet_nupfreq) * uint64_t(packet_ntsamp);
	if (_unlikely(!n4 || (n4 > 10000)))
	    continue;    // silently drop bad packets
	if (_unlikely(data_nbytes != int(n4)))
	    continue;    // silently drop bad packets
	if (_unlikely(packet_nbytes != packet_size(packet_nbeam, packet_nfreq, packet_nupfreq, packet_ntsamp)))
	    continue;    // silently drop bad packets
	
	const uint16_t *packet_beam_ids = (const uint16_t *) (packet + 24);
	const uint8_t *packet_freq_ids = packet + 24 + 2*packet_nbeam;
	const uint8_t *packet_scales = packet + 24 + 2*packet_nbeam + 2*packet_nfreq;
	const uint8_t *packet_offsets = packet + 24 + 2*packet_nbeam + 2*packet_nfreq + 4*packet_nbeam*packet_nfreq;
	const uint8_t *packet_data = packet + 24 + 2*packet_nbeam + 2*packet_nfreq + 8*packet_nbeam*packet_nfreq;

	if (npackets_received == 0) {
	    for (int i = 0; i < nassemblers; i++)
		assemblers[i]->start_stream(packet_fpga_counts_per_sample, packet_nupfreq);
	}

	npackets_received++;

	// Loop over beams in the packet, matching to beam_assembler objects.

	for (int ibeam = 0; ibeam < packet_nbeam; ibeam++) {

	    // Find assembler matching beam ID	    
	    int beam_id = packet_beam_ids[ibeam];

	    int assembler_ix;
	    for (assembler_ix = 0; assembler_ix < nassemblers; assembler_ix++) {
		if (assembler_beam_ids[assembler_ix] == beam_id)
		    break;   // found assembler
	    }

	    if (assembler_ix >= nassemblers)
		continue;   // no assembler found, proceed to next beam_id
	    
	    // Make subpacket
	    
	    uint8_t *subpacket = assembler_packet_lists[assembler_ix].data_end;
	    uint8_t *subpacket_freq_ids = subpacket + 26;
	    uint8_t *subpacket_scales = subpacket + 26 + 2*packet_nfreq;
	    uint8_t *subpacket_offsets = subpacket + 26 + 6*packet_nfreq;
	    uint8_t *subpacket_data = subpacket + 26 + 10*packet_nfreq;
	    
	    int subpacket_data_size = packet_nfreq * packet_nupfreq * packet_ntsamp;
	    int subpacket_total_size = subpacket_data_size + 26 + 10*packet_nfreq;
	    
	    memcpy(subpacket, packet, 24);                            // copy header and overwrite a few fields...
	    *((int16_t *) (subpacket+4)) = subpacket_data_size;       // overwrite 'data_nbytes' field
	    *((uint16_t *) (subpacket+16)) = uint16_t(1);             // overwrite 'nbeams' field
	    *((uint16_t *) (subpacket+24)) = packet_beam_ids[ibeam];  // beam_id field
	    
	    memcpy(subpacket_freq_ids, packet_freq_ids, 2*packet_nfreq);
	    memcpy(subpacket_scales, packet_scales + 4*ibeam*packet_nfreq, 4*packet_nfreq);
	    memcpy(subpacket_offsets, packet_offsets + 4*ibeam*packet_nfreq, 4*packet_nfreq);
	    memcpy(subpacket_data, packet_data + ibeam*subpacket_data_size, subpacket_data_size);

	    // assembler_timestamp gets set when first packet is added
	    if (assembler_packet_lists[assembler_ix].curr_npackets == 0)
		assembler_timestamps[assembler_ix] = curr_timestamp;

	    assembler_packet_lists[assembler_ix].add_packet(subpacket_total_size);

	    if (!assembler_packet_lists[assembler_ix].is_full)
		continue;

	    // Packet list is full, so send it to the assembler thread.
	    bool assembler_alive = assemblers[assembler_ix]->put_unassembled_packets(assembler_packet_lists[assembler_ix]);

	    if (!assembler_alive) {
		// Is this what we should do?
		cerr << "ch_frb_io:  assembler thread died unexpectedly!  network thread will die too...\n";
		return npackets_received;
	    }
	}
    }

    for (int i = 0; i < nassemblers; i++)
	if (assembler_packet_lists[i].curr_npackets > 0)
	    assemblers[i]->put_unassembled_packets(assembler_packet_lists[i]);

    return npackets_received;
}


}  // namespace ch_frb_io
