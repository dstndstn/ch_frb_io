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
static ssize_t network_thread_main2(intensity_network_stream *stream, udp_packet_list *packet_lists);   // returns npackets_received


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


void intensity_network_stream::end_stream(bool join_threads)
{
    bool call_join_after_releasing_lock = false;

    pthread_mutex_lock(&this->lock);

    // Set flags as if stream had run to completion.  This is convenient e.g. for waking up threads
    // which are blocked in wait_for_packets().
    this->stream_started = true;
    this->stream_ended = true;

    if (join_threads && !network_thread_joined) {
	call_join_after_releasing_lock = true;
	this->network_thread_joined = true;
    }
    
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->lock);

    for (unsigned int i = 0; i < assemblers.size(); i++)
	assemblers[i]->end_stream(join_threads);

    if (call_join_after_releasing_lock)
	pthread_join(network_thread, NULL);
}


void intensity_network_stream::wait_for_end_of_stream(bool join_threads)
{
    pthread_mutex_lock(&lock);
    
    if (!stream_started) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::wait_for_end_of_stream() was called with no prior call to start_stream()");
    }
    
    while (!stream_ended)
	pthread_cond_wait(&this->cond_state_changed, &this->lock);

    pthread_mutex_unlock(&lock);

    // Looks weird but does the correct thing
    this->end_stream(join_threads);
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


// FIXME also in intensity_network_ostream.cpp
inline int packet_size(int nbeam, int nfreq, int nupfreq, int ntsamp)
{
    return 24 + 2*nbeam + 2*nfreq + 8*nbeam*nfreq + (nbeam * nfreq * nupfreq * ntsamp);
}


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

    cerr << "ch_frb_io: network thread starting\n";
    stream->network_thread_startup();

    ssize_t npackets_received = 0;

    try {
	// It's convenient to allocate the udp_packet_lists before calling network_thread_main2()
	// FIXME this seems a little awkward, is there a better way?

	int nassemblers = stream->assemblers.size();

	udp_packet_list *packet_lists = new udp_packet_list[nassemblers];
	for (int i = 0; i < nassemblers; i++)
	    packet_lists[i].initialize(stream->assemblers[i]->beam_id);
	
	npackets_received = network_thread_main2(stream.get(), packet_lists);

	for (int i = 0; i < nassemblers; i++)
	    packet_lists[i].destroy();

	delete[] packet_lists;

    } catch (...) {
	stream->end_stream(false);   // "false" means "don't join threads" (would deadlock otherwise!)
	throw;
    }

    stream->end_stream(false);   // "false" has same meaning as above

    cerr << ("ch_frb_io: network thread exiting (" + to_string(npackets_received) + " packets received)\n");
    return NULL;
}


// Returns number of packets received
static ssize_t network_thread_main2(intensity_network_stream *stream, udp_packet_list *packet_lists)
{
    // FIXME is 2MB socket_bufsize a good choice?  I would have guessed a larger value 
    // would be better, but 2MB is the max allowed on my osx laptop.
    static const int socket_bufsize = 2 << 21; 

    int nassemblers = stream->assemblers.size();

    intensity_beam_assembler *assemblers[nassemblers];
    for (int i = 0; i < nassemblers; i++)
	assemblers[i] = stream->assemblers[i].get();

    vector<uint8_t> packet_vec(udp_packet_list::max_packet_size);
    uint8_t *packet = &packet_vec[0];

    //
    // Create socket
    //

    int sock_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (sock_fd < 0)
	throw runtime_error(string("ch_frb_io: socket() failed: ") + strerror(errno));

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));
	
    server_address.sin_family = AF_INET;
    inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
    server_address.sin_port = htons(stream->udp_port);
    
    int err = ::bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: bind() failed: ") + strerror(errno));

    err = setsockopt(sock_fd, SOL_SOCKET, SO_RCVBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt() failed: ") + strerror(errno));

    cerr << ("ch_frb_io: network thread: listening for packets on port " + to_string(stream->udp_port) + "\n");

    //
    // Wait for stream to start
    //

    bool cancelled = !stream->wait_for_start_stream();
    
    if (cancelled)
	return 0;
    
    //
    // Main packet loop!
    // 
    // FIXME there is a lot of silent "swallowing" of error conditions here, needs revisiting to
    // think about what really makes sense.
    //

    ssize_t npackets_received = 0;

    for (;;) {
	int packet_nbytes = read(sock_fd, (char *) packet, udp_packet_list::max_packet_size);

	if (_unlikely(packet_nbytes < 0))
	    throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
	if (_unlikely(packet_nbytes >= udp_packet_list::max_packet_size))
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
	    return npackets_received;   // Note: this branch is the only way out of the main packet loop
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
		if (packet_lists[assembler_ix].beam_id == beam_id)
		    break;   // found assembler
	    }

	    if (assembler_ix >= nassemblers)
		continue;   // no assembler found, proceed to next beam_id
	    
	    // Make subpacket
	    
	    uint8_t *subpacket = packet_lists[assembler_ix].data_end;
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
	    
	    packet_lists[assembler_ix].add_packet(subpacket_total_size);

	    if (!packet_lists[assembler_ix].is_full)
		continue;

	    // Packet list is full, so send it to the assembler thread.
	    bool assembler_alive = assemblers[assembler_ix]->put_unassembled_packets(packet_lists[assembler_ix]);

	    if (!assembler_alive) {
		// Is this what we should do?
		cerr << "ch_frb_io:  assembler thread died unexpectedly!  network thread will die too...\n";
		return npackets_received;
	    }
	}
    }
}


}  // namespace ch_frb_io
