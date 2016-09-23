#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// class intensity_network_stream


// Static member function (de facto constructor)
shared_ptr<intensity_network_stream> intensity_network_stream::make(const initializer &x)
{
    intensity_network_stream *retp = new intensity_network_stream(x);
    shared_ptr<intensity_network_stream> ret(retp);

    ret->_open_socket();

    // Spawn assembler thread
    int err = pthread_create(&ret->assembler_thread, NULL, assembler_pthread_main, (void *) &ret);
    if (err)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_stream constructor: ") + strerror(errno));
    
    // Wait for assembler thread to start
    pthread_mutex_lock(&ret->state_lock);
    while (!ret->assembler_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->state_lock);
    pthread_mutex_unlock(&ret->state_lock);    

    // Spawn network thread
    err = pthread_create(&ret->network_thread, NULL, network_pthread_main, (void *) &ret);
    if (err)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_stream constructor: ") + strerror(errno));
    
    // Wait for network thread to start
    pthread_mutex_lock(&ret->state_lock);
    while (!ret->network_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->state_lock);
    pthread_mutex_unlock(&ret->state_lock);    

    return ret;
}


intensity_network_stream::intensity_network_stream(const initializer &x) :
    ini_params(x),
    nassemblers(x.beam_ids.size())
{
    // Argument checking

    if (nassemblers == 0)
	throw runtime_error("ch_frb_io: length-zero beam_id vector passed to intensity_network_stream constructor");

    for (int i = 0; i < nassemblers; i++) {
	if ((x.beam_ids[i] < 0) || (x.beam_ids[i] > constants::max_allowed_beam_id))
	    throw runtime_error("ch_frb_io: bad beam_id passed to intensity_network_stream constructor");
	for (int j = 0; j < i; j++)
	    if (x.beam_ids[i] == x.beam_ids[j])
		throw runtime_error("ch_frb_io: duplicate beam_ids passed to intensity_network_stream constructor");
    }

    if ((x.udp_port <= 0) || (x.udp_port >= 65536))
	throw runtime_error("ch_frb_io: intensity_network_stream constructor: bad udp port " + to_string(x.udp_port));

    if (x.mandate_fast_kernels && x.mandate_reference_kernels)
	throw runtime_error("ch_frb_io: both flags mandate_fast_kernels, mandate_reference_kernels were set");

    // All initializations except the socket (which is initialized in _open_socket()),
    // and the assemblers (which are initalized when the first packet is received).

    int capacity = constants::unassembled_ringbuf_capacity;
    int max_npackets = constants::max_unassembled_packets_per_list;
    int max_nbytes = constants::max_unassembled_nbytes_per_list;
    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (capacity, max_npackets, max_nbytes);

    this->incoming_packet_list = udp_packet_list(constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list);

    this->event_counts = vector<int64_t> (event_type::num_types, 0);
    this->network_thread_event_subcounts = vector<int64_t> (event_type::num_types, 0);
    this->assembler_thread_event_subcounts = vector<int64_t> (event_type::num_types, 0);

    pthread_mutex_init(&state_lock, NULL);
    pthread_mutex_init(&event_lock, NULL);
    pthread_cond_init(&cond_state_changed, NULL);
}


intensity_network_stream::~intensity_network_stream()
{
    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&state_lock);
    pthread_mutex_destroy(&event_lock);

    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }
}


// Socket initialization factored to its own routine, rather than putting it in the constructor,
// so that the socket will always be closed if an exception is thrown somewhere.
void intensity_network_stream::_open_socket()
{
    const int socket_bufsize = constants::recv_socket_bufsize;
    const struct timeval tv_timeout = { 0, constants::recv_socket_timeout_usec };

    this->sockfd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (sockfd < 0)
	throw runtime_error(string("ch_frb_io: socket() failed: ") + strerror(errno));

    int err = setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVBUF) failed: ") + strerror(errno));

    err = setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, &tv_timeout, sizeof(tv_timeout));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVTIMEO) failed: ") + strerror(errno));
}


void intensity_network_stream::_add_event_counts(vector<int64_t> &event_subcounts)
{
    if (event_counts.size() != event_subcounts.size())
	throw runtime_error("ch_frb_io: internal error: vector length mismatch in intensity_network_stream::_add_event_counts()");

    pthread_mutex_lock(&this->event_lock);
    for (unsigned int i = 0; i < event_counts.size(); i++)
	event_counts[i] += event_subcounts[i];
    pthread_mutex_unlock(&this->event_lock);

    memset(&event_subcounts[0], 0, event_subcounts.size() * sizeof(event_subcounts[0]));
}


void intensity_network_stream::start_stream()
{
    pthread_mutex_lock(&this->state_lock);

    if (stream_started) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::start_stream() called on running, completed, or cancelled stream");
    }

    this->stream_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);
}


// Just sets the stream_ended flag and returns.
// The network thread will see that this flag has been set, flush packets to the assembler thread, call unassembled_ringbuf->end_stream(), and exit.
// The assembler thread will see that the ringbuf has ended, flush buffers to the processing threads, and exit.
void intensity_network_stream::end_stream()
{
    pthread_mutex_lock(&this->state_lock);
    this->stream_started = true;
    this->first_packet_received = true;
    this->stream_ended = true;    
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);
}


void intensity_network_stream::join_threads()
{
    pthread_mutex_lock(&this->state_lock);
    
    if (!stream_started) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::join_network_thread() was called with no prior call to start_stream()");
    }

    if (join_called) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: double call to intensity_network_stream::join_network_thread()");
    }

    this->join_called = true;
    pthread_mutex_unlock(&this->state_lock);

    pthread_join(network_thread, NULL);
    pthread_join(assembler_thread, NULL);
}


shared_ptr<assembled_chunk> intensity_network_stream::get_assembled_chunk(int assembler_index)
{
    if ((assembler_index < 0) || (assembler_index >= nassemblers))
	throw runtime_error("ch_frb_io: bad assembler_ix passed to intensity_network_stream::get_assembled_chunk()");

    // Wait for first_packet_received flag to be set.
    pthread_mutex_lock(&this->state_lock);
    while (!first_packet_received)
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    pthread_mutex_unlock(&this->state_lock);

    // There is a corner case where the vector is still length-zero after the flag gets set.
    // This happens if the stream was asynchronous cancelled before receiving the first packet.
    if (assemblers.size() == 0)
	return shared_ptr<assembled_chunk> ();

    return assemblers[assembler_index]->get_assembled_chunk();
}


intensity_network_stream::initializer intensity_network_stream::get_initializer()
{
    return this->ini_params;
}


vector<int64_t> intensity_network_stream::get_event_counts()
{
    vector<int64_t> ret(event_type::num_types, 0);

    pthread_mutex_lock(&this->event_lock);
    memcpy(&ret[0], &event_counts[0], ret.size() * sizeof(ret[0]));
    pthread_mutex_unlock(&this->event_lock);    

    return ret;
}


bool intensity_network_stream::get_first_packet_params(int &nupfreq, int &nt_per_packet, uint64_t &fpga_counts_per_sample, uint64_t &fpga_count)
{
    pthread_mutex_lock(&this->state_lock);

    while (!this->first_packet_received)
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    
    if (this->stream_ended) {
	pthread_mutex_unlock(&this->state_lock);
	nupfreq = 0;
	nt_per_packet = 0;
	fpga_counts_per_sample = 0;
	fpga_count = 0;
	return false;
    }
    
    pthread_mutex_unlock(&this->state_lock);

    // OK to make these assignments without holding lock
    nupfreq = fp_nupfreq;
    nt_per_packet = fp_nt_per_packet;
    fpga_counts_per_sample = fp_fpga_counts_per_sample;
    fpga_count = fp_fpga_count;
    return true;
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


// static member function
void *intensity_network_stream::network_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_pthread_main()");

    shared_ptr<intensity_network_stream> *arg = (shared_ptr<intensity_network_stream> *) opaque_arg;
    shared_ptr<intensity_network_stream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to network_pthread_main()");

    stream->_network_thread_start();

    // We use try..catch to ensure that _network_thread_exit() always gets called, even if an exception is thrown.
    try {
	stream->_network_thread_body();
    } catch (...) {
	stream->_network_thread_exit();
	throw;
    }
    
    stream->_network_thread_exit();
    return NULL;
}


void intensity_network_stream::_network_thread_start()
{
    // Advance stream state to "network_thread_started" state,
    // and wait for another thread to advance it to "stream_started" state.

    pthread_mutex_lock(&this->state_lock);

    this->network_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);

    for (;;) {
	if (this->stream_ended) {
	    pthread_mutex_unlock(&this->state_lock);
	    return;
	}
	if (this->stream_started) {
	    pthread_mutex_unlock(&this->state_lock);
	    break;
	}
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    }
}


void intensity_network_stream::_network_thread_body()
{
    // Start listening on socket 

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));
	
    server_address.sin_family = AF_INET;
    inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
    server_address.sin_port = htons(ini_params.udp_port);

    int err = ::bind(sockfd, (struct sockaddr *) &server_address, sizeof(server_address));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: bind() failed: ") + strerror(errno));

    cerr << ("ch_frb_io: listening for packets on port " + to_string(ini_params.udp_port) + "\n");

    // Main packet loop

    int64_t *event_subcounts = &this->network_thread_event_subcounts[0];
    struct timeval tv_ini = xgettimeofday();
    uint64_t incoming_packet_list_timestamp = 0;
    uint64_t cancellation_check_timestamp = 0;
    bool is_first_packet = true;

    for (;;) {
	// All timestamps are in microseconds relative to tv_ini.
	uint64_t curr_timestamp = usec_between(tv_ini, xgettimeofday());

	// Periodically check whether stream has been cancelled by end_stream().
	if (curr_timestamp > cancellation_check_timestamp + constants::stream_cancellation_latency_usec) {
	    pthread_mutex_lock(&this->state_lock);

	    if (this->stream_ended) {
		pthread_mutex_unlock(&this->state_lock);    
		return;
	    }

	    pthread_mutex_unlock(&this->state_lock);
	    cancellation_check_timestamp = curr_timestamp;
	}

	// Periodically flush packets to assembler thread (only happens if packet rate is low; normal case is that the packet_list fills first)
	if (curr_timestamp > incoming_packet_list_timestamp + constants::unassembled_ringbuf_timeout_usec)
	    this->_put_unassembled_packets();   // Note: calls _add_event_counts()

	// Read new packet from socket (note that socket has a timeout, so this call can time out)
	uint8_t *packet_data = incoming_packet_list.data_end;
	int packet_nbytes = read(sockfd, packet_data, constants::max_input_udp_packet_size + 1);

	// Check for error or timeout in read()
	if (packet_nbytes < 0) {
	    if ((errno == EAGAIN) || (errno == ETIMEDOUT))
		continue;  // normal timeout
	    throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
	}

	event_subcounts[event_type::byte_received] += packet_nbytes;
	event_subcounts[event_type::packet_received]++;

	// If we receive a special "short" packet (length 24), it indicates end-of-stream.
	// FIXME is this a temporary kludge or something which should be documented in the packet protocol?
	if (_unlikely(packet_nbytes == 24)) {
	    event_subcounts[event_type::packet_end_of_stream]++;
	    if (ini_params.accept_end_of_stream_packets)
		return;   // triggers shutdown of entire stream
	    continue;
	}

	// Special logic for processing first packet, which triggers some initializations.
	if (is_first_packet) {
	    intensity_packet packet;

	    if (!packet.read(packet_data, packet_nbytes)) {
		event_subcounts[event_type::packet_bad]++;
		continue;
	    }

	    this->fp_nupfreq = packet.nupfreq;
	    this->fp_nt_per_packet = packet.ntsamp;
	    this->fp_fpga_counts_per_sample = packet.fpga_counts_per_sample;
	    this->fp_fpga_count = packet.fpga_count;
	    is_first_packet = false;

	    // Now that we know the first packet params, we can initialize the assemblers
	    this->assemblers.resize(nassemblers);
	    for (int ix = 0; ix < nassemblers; ix++)
		assemblers[ix] = make_unique<assembled_chunk_ringbuf>(ini_params, ini_params.beam_ids[ix], fp_nupfreq, 
								      fp_nt_per_packet, fp_fpga_counts_per_sample, fp_fpga_count);

	    // This will unblock the assembler thread, which waits for the first_packet_received
	    // flag as soon as it is spawned.
	    pthread_mutex_lock(&this->state_lock);
	    this->first_packet_received = true;
	    pthread_cond_broadcast(&this->cond_state_changed);
	    pthread_mutex_unlock(&this->state_lock);
	}

	// The incoming_packet_list is timestamped with the arrival time of its first packet.
	if (incoming_packet_list.curr_npackets == 0)
	    incoming_packet_list_timestamp = curr_timestamp;

	incoming_packet_list.add_packet(packet_nbytes);

	if (incoming_packet_list.is_full)
	    this->_put_unassembled_packets();
    }
}


// This gets called when the network thread exits, either 
// because the stream ended or an exception was thrown.
void intensity_network_stream::_network_thread_exit()
{
    // This just sets the stream_ended flag, if it hasn't been set yet.
    this->end_stream();

    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }
    
    // Note: calls _add_event_counts()
    this->_put_unassembled_packets();

    unassembled_ringbuf->end_stream();
}


void intensity_network_stream::_put_unassembled_packets()
{
    int npackets = incoming_packet_list.curr_npackets;
    if (!npackets)
	return;

    bool success = unassembled_ringbuf->put_packet_list(incoming_packet_list, false);

    if (!success) {
	network_thread_event_subcounts[event_type::packet_dropped] += npackets;

	if (ini_params.emit_warning_on_buffer_drop)
	    cerr << "ch_frb_io: assembler thread crashed or is running slow, dropping packets\n";
	if (ini_params.throw_exception_on_buffer_drop)
	    throw runtime_error("ch_frb_io: packets were dropped and stream was constructed with 'throw_exception_on_buffer_drop' flag");
    }

    this->_add_event_counts(network_thread_event_subcounts);
}


// -------------------------------------------------------------------------------------------------
//
// assembler thread


// static member function
void *intensity_network_stream::assembler_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to assembler_pthread_main()");

    shared_ptr<intensity_network_stream> *arg = (shared_ptr<intensity_network_stream> *) opaque_arg;
    shared_ptr<intensity_network_stream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to assembler_thread_main()");

    stream->_assembler_thread_start();

    // We use try..catch to ensure that _assembler_thread_exit() always gets called, even if an exception is thrown.
    try {
	stream->_assembler_thread_body();
    } catch (...) {
	stream->_assembler_thread_exit();
	throw;
    }

    stream->_assembler_thread_exit();
    return NULL;
}


void intensity_network_stream::_assembler_thread_start()
{
    // Set the assembler_thread_started flag (this unblocks the thread which spawned the assembler thread)
    pthread_mutex_lock(&this->state_lock);
    this->assembler_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);

    // Wait for first_packet_received flag to be set.
    pthread_mutex_lock(&this->state_lock);
    while (!first_packet_received)
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    pthread_mutex_unlock(&this->state_lock);
}


void intensity_network_stream::_assembler_thread_body()
{
    udp_packet_list packet_list(constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list);
    int64_t *event_subcounts = &this->assembler_thread_event_subcounts[0];

    while (unassembled_ringbuf->get_packet_list(packet_list)) {
	for (int ipacket = 0; ipacket < packet_list.curr_npackets; ipacket++) {
            uint8_t *packet_data = packet_list.get_packet_data(ipacket);
            int packet_nbytes = packet_list.get_packet_nbytes(ipacket);
	    intensity_packet packet;

	    if (!packet.read(packet_data, packet_nbytes)) {
		event_subcounts[event_type::packet_bad]++;
		continue;
	    }

	    bool mismatch = (packet.nupfreq != fp_nupfreq) || (packet.ntsamp != fp_nt_per_packet) || (packet.fpga_counts_per_sample != fp_fpga_counts_per_sample);

	    if (_unlikely(mismatch)) {
		event_subcounts[event_type::first_packet_mismatch]++;
		continue;
	    }

	    // All checks passed.  Packet is declared "good" here.  
	    //
	    // The following checks have been performed, either in this routine or in intensity_packet::read().
	    //   - dimensions (nbeams, nfreq_coarse, nupfreq, ntsamp) are positive,
	    //     and not large enough to lead to integer overflows
	    //   - packet and data byte counts are correct
	    //   - coarse_freq_ids are valid (didn't check for duplicates but that's ok)
	    //   - ntsamp is a power of two
	    //   - fpga_counts_per_sample is > 0
	    //   - fpga_count is a multiple of (fpga_counts_per_sample * ntsamp)
	    //
	    // These checks are assumed by assembled_chunk::add_packet(), and mostly aren't rechecked, 
	    // so it's important that they're done here!

	    event_subcounts[event_type::packet_good]++;
	    
	    int nbeams = packet.nbeams;
	    int nfreq_coarse = packet.nfreq_coarse;
	    int new_data_nbytes = nfreq_coarse * packet.nupfreq * packet.ntsamp;
	    const int *assembler_beam_ids = &ini_params.beam_ids[0];  // bare pointer for speed
    
	    // Danger zone: we modify the packet by leaving its pointers in place, but shortening its
	    // length fields.  The new packet corresponds to a subset of the original packet containing
	    // only beam index zero.
	    
	    packet.data_nbytes = new_data_nbytes;
	    packet.nbeams = 1;
	    
	    for (int ibeam = 0; ibeam < nbeams; ibeam++) {
		// Loop invariant: at the top of this loop, 'packet' corresponds to a subset of the
		// original packet containing only beam index 'ibeam'.
		
		// Loop over assembler ids, to find a match for the packet_id.
		int packet_id = packet.beam_ids[0];
		int assembler_ix = 0;

		for (;;) {
		    if (assembler_ix >= nassemblers) {
			// No match found
			event_subcounts[event_type::beam_id_mismatch]++;
			break;
		    }

		    if (assembler_beam_ids[assembler_ix] != packet_id) {
			assembler_ix++;
			continue;
		    }

		    // Match found
		    assemblers[assembler_ix]->put_unassembled_packet(packet, event_subcounts);
		    break;
		}
		
		// Danger zone: we do some pointer arithmetic, to modify the packet so that it now
		// corresponds to a new subset of the original packet, corresponding to beam index (ibeam+1).
		
		packet.beam_ids += 1;
		packet.scales += nfreq_coarse;
		packet.offsets += nfreq_coarse;
		packet.data += new_data_nbytes;
	    }
	}

	this->_add_event_counts(assembler_thread_event_subcounts);
    }
}


void intensity_network_stream::_assembler_thread_exit()
{
    unassembled_ringbuf->end_stream();

    // Use assemblers.size() instead of 'nassemblers', to correctly handle a corner case where the stream
    // is asynchronously cancelled before receiving the first packet, and the vector never gets allocated.

    for (unsigned int i = 0; i < assemblers.size(); i++)
	assemblers[i]->end_stream(&assembler_thread_event_subcounts[0]);

    this->_add_event_counts(assembler_thread_event_subcounts);

    vector<int64_t> counts = this->get_event_counts();

    stringstream ss;
    ss << "ch_frb_io: assembler thread exiting\n"
       << "    bytes recevied (GB): " << (1.0e-9 * counts[event_type::byte_received]) << "\n"
       << "    packets received: " << counts[event_type::packet_received] << "\n"
       << "    good packets: " << counts[event_type::packet_good] << "\n"
       << "    bad packets: " << counts[event_type::packet_bad] << "\n"
       << "    dropped packets: " << counts[event_type::packet_dropped] << "\n"
       << "    end-of-stream packets: " << counts[event_type::packet_end_of_stream] << "\n"
       << "    beam id mismatches: " << counts[event_type::beam_id_mismatch] << "\n"
       << "    first-packet mismatches: " << counts[event_type::first_packet_mismatch] << "\n"
       << "    assembler hits: " << counts[event_type::assembler_hit] << "\n"
       << "    assembler misses: " << counts[event_type::assembler_miss] << "\n"
       << "    assembled chunks dropped: " << counts[event_type::assembled_chunk_dropped] << "\n"
       << "    assembled chunks queued: " << counts[event_type::assembled_chunk_queued] << "\n";

    cerr << ss.str().c_str();
}


}  // namespace ch_frb_io
