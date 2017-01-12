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

    // Spawn assembler thread.  Note that we pass a bare pointer to an object ('ret') on our stack
    // and this pointer will be invalid after make() returns.  Therefore, the assembler thread only
    // dereferences the pointer before it sets the assembler_thread_started flag, and make() waits for this
    // flag to be set before it returns.  The same synchronization logic applies to the network thread.

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


intensity_network_stream::intensity_network_stream(const initializer &ini_params_) :
    ini_params(ini_params_),
    nassemblers(ini_params_.beam_ids.size())
{
    // Argument checking

    if (nassemblers == 0)
	throw runtime_error("ch_frb_io: length-zero beam_id vector passed to intensity_network_stream constructor");

    for (int i = 0; i < nassemblers; i++) {
	if ((ini_params.beam_ids[i] < 0) || (ini_params.beam_ids[i] > constants::max_allowed_beam_id))
	    throw runtime_error("ch_frb_io: bad beam_id passed to intensity_network_stream constructor");
	for (int j = 0; j < i; j++)
	    if (ini_params.beam_ids[i] == ini_params.beam_ids[j])
		throw runtime_error("ch_frb_io: duplicate beam_ids passed to intensity_network_stream constructor");
    }

    if ((ini_params.udp_port <= 0) || (ini_params.udp_port >= 65536))
	throw runtime_error("ch_frb_io: intensity_network_stream constructor: bad udp port " + to_string(ini_params.udp_port));

    if (ini_params.mandate_fast_kernels && ini_params.mandate_reference_kernels)
	throw runtime_error("ch_frb_io: both flags mandate_fast_kernels, mandate_reference_kernels were set");

#ifndef __AVX2__
    if (ini_params.mandate_fast_kernels)
	throw runtime_error("ch_frb_io: the 'mandate_fast_kernels' flag was set, but this machine does not have the AVX2 instruction set");
#endif

    // All initializations except the socket (which is initialized in _open_socket()),
    // and the assemblers (which are initalized when the first packet is received).

    int capacity = constants::unassembled_ringbuf_capacity;
    int max_npackets = constants::max_unassembled_packets_per_list;
    int max_nbytes = constants::max_unassembled_nbytes_per_list;
    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (capacity, max_npackets, max_nbytes);

    this->incoming_packet_list = make_unique<udp_packet_list> (constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list);

    this->cumulative_event_counts = vector<int64_t> (event_type::num_types, 0);
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

    // bufsize
    int err = setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVBUF) failed: ") + strerror(errno));

    // timeout
    err = setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, &tv_timeout, sizeof(tv_timeout));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVTIMEO) failed: ") + strerror(errno));
}


void intensity_network_stream::_add_event_counts(vector<int64_t> &event_subcounts)
{
    if (cumulative_event_counts.size() != event_subcounts.size())
	throw runtime_error("ch_frb_io: internal error: vector length mismatch in intensity_network_stream::_add_event_counts()");

    pthread_mutex_lock(&this->event_lock);
    for (unsigned int i = 0; i < cumulative_event_counts.size(); i++)
	this->cumulative_event_counts[i] += event_subcounts[i];
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


// Just sets the stream_end_requested flag and returns.  The shutdown logic then proceeeds as follows.
// The network thread will see that the stream_end_requested flag has been set, flush packets to the 
// assembler thread, call unassembled_ringbuf->end_stream(), and exit.  The assembler thread will see 
// that the ringbuf has ended, flush assembled_chunks to the processing threads, and exit.
//
// Note that end_stream() can be called multiple times (this usually happens as part of the shutdown process).

void intensity_network_stream::end_stream()
{
    pthread_mutex_lock(&this->state_lock);
    this->stream_started = true;
    this->first_packet_received = true;
    this->assemblers_initialized = true;
    this->stream_end_requested = true;    
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);
}


void intensity_network_stream::join_threads()
{
    pthread_mutex_lock(&this->state_lock);
    
    if (!stream_started) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::join_threads() was called with no prior call to start_stream()");
    }

    if (join_called) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: double call to intensity_network_stream::join_threads()");
    }

    this->join_called = true;
    pthread_mutex_unlock(&this->state_lock);

    int err = pthread_join(network_thread, NULL);
    if (err)
	throw runtime_error("ch_frb_io: couldn't join network thread [input]");

    err = pthread_join(assembler_thread, NULL);
    if (err)
	throw runtime_error("ch_frb_io: couldn't join assembler thread");
}


shared_ptr<assembled_chunk> intensity_network_stream::get_assembled_chunk(int assembler_index)
{
    if ((assembler_index < 0) || (assembler_index >= nassemblers))
	throw runtime_error("ch_frb_io: bad assembler_ix passed to intensity_network_stream::get_assembled_chunk()");

    // Wait for assemblers_initialized flag to be set.
    // After this, it's OK to access 'assemblers' without holding the lock.
    pthread_mutex_lock(&this->state_lock);
    while (!this->assemblers_initialized)
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    pthread_mutex_unlock(&this->state_lock);

    // This can happen if end_stream() was called before the assemblers were initialized.
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
    memcpy(&ret[0], &this->cumulative_event_counts[0], ret.size() * sizeof(ret[0]));
    pthread_mutex_unlock(&this->event_lock);    

    return ret;
}

unordered_map<string, uint64_t> intensity_network_stream::get_perhost_packets()
{
    pthread_mutex_lock(&this->event_lock);
    unordered_map<uint64_t, uint64_t> raw(perhost_packets);
    pthread_mutex_unlock(&this->event_lock);
    // Convert to strings
    unordered_map<string, uint64_t> rtn;
    for (auto it = raw.begin(); it != raw.end(); it++) {
        // IPv4 address in high 32 bits, port in low 16 bits
        uint32_t ip = (it->first >> 32) & 0xffffffff;
        uint32_t port = (it->first & 0xffff);
        string sender = to_string(ip >> 24) + "." + to_string((ip >> 16) & 0xff)
            + "." + to_string((ip >> 8) & 0xff) + "." + to_string(ip & 0xff)
            + ":" + to_string(port);
        rtn[sender] = it->second;
    }
    return rtn;
}

vector<unordered_map<string, uint64_t> >
intensity_network_stream::get_statistics() {
    vector<unordered_map<string, uint64_t> > R;
    bool first_packet;
    bool assemblers_init;
    unordered_map<string, uint64_t> m;

    pthread_mutex_lock(&this->state_lock);
    first_packet = this->first_packet_received;
    assemblers_init = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);

    // Collect statistics for this stream as a whole:
    m["first_packet_received"]  =  first_packet;
    m["nupfreq"]                = (first_packet ? this->fp_nupfreq : 0);
    m["nt_per_packet"]          = (first_packet ? this->fp_nt_per_packet : 0);
    m["fpga_counts_per_sample"] = (first_packet ? this->fp_fpga_counts_per_sample : 0);
    m["fpga_count"]             = (first_packet ? this->fp_fpga_count : 0);

    vector<int64_t> counts = get_event_counts();
    m["count_bytes_received"     ] = counts[event_type::byte_received];
    m["count_packets_received"   ] = counts[event_type::packet_received];
    m["count_packets_good"       ] = counts[event_type::packet_good];
    m["count_packets_bad"        ] = counts[event_type::packet_bad];
    m["count_packets_dropped"    ] = counts[event_type::packet_dropped];
    m["count_packets_endofstream"] = counts[event_type::packet_end_of_stream];
    m["count_beam_id_mismatch"   ] = counts[event_type::beam_id_mismatch];
    m["count_stream_mismatch"    ] = counts[event_type::first_packet_mismatch];
    m["count_assembler_hits"     ] = counts[event_type::assembler_hit];
    m["count_assembler_misses"   ] = counts[event_type::assembler_miss];
    m["count_assembler_drops"    ] = counts[event_type::assembled_chunk_dropped];
    m["count_assembler_queued"   ] = counts[event_type::assembled_chunk_queued];

    int nbeams = this->ini_params.beam_ids.size();
    m["nbeams"] = nbeams;
    R.push_back(m);

    // Report per-host packet counts
    R.push_back(this->get_perhost_packets());

    // Collect statistics per beam:
    for (int b=0; b<nbeams; b++) {
        m.clear();
        m["beam_id"] = this->ini_params.beam_ids[b];
        if (assemblers_init && (assemblers.size() > 0)) {
            // Grab the ring buffer to find the min & max chunk numbers and size.
            uint64_t ichunk, nchunks, capacity, ntot, oldest;
            this->assemblers[b]->get_ringbuf_size(&ichunk, &nchunks, &capacity, &ntot, &oldest);
            m["ringbuf_next"] = ichunk;
            m["ringbuf_ready"] = nchunks;
            m["ringbuf_capacity"] = capacity;
            m["ringbuf_ntotal"] = ntot;
            if (ntot == 0) {
                m["ringbuf_chunk_min"] = 0;
                m["ringbuf_chunk_max"] = 0;
            } else {
                m["ringbuf_chunk_min"] = oldest;
                m["ringbuf_chunk_max"] = ichunk + nchunks - 1;
            }
        }
        R.push_back(m);
    }

    return R;
}

vector< vector< shared_ptr<assembled_chunk> > >
intensity_network_stream::get_ringbuf_snapshots(vector<uint64_t> beams)
{
    vector< vector< shared_ptr<assembled_chunk> > > R;

    pthread_mutex_lock(&this->state_lock);
    bool assemblers_init = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);

    if (!assemblers_init)
        return R;

    int nbeams = this->ini_params.beam_ids.size();

    for (int ib=0; ib<beams.size(); ib++) {
        uint64_t beam = beams[ib];
        vector<shared_ptr<assembled_chunk> > chunks;
        // Which of my assemblers (if any) is handling the requested beam?
        for (int i=0; i<nbeams; i++) {
            if (this->ini_params.beam_ids[i] != beam)
                continue;
            chunks = this->assemblers[i]->get_ringbuf_snapshot();
        }
        R.push_back(chunks);
    }
    return R;
}

bool intensity_network_stream::get_first_packet_params(int &nupfreq, int &nt_per_packet, uint64_t &fpga_counts_per_sample, uint64_t &fpga_count)
{
    pthread_mutex_lock(&this->state_lock);

    while (!this->first_packet_received)
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    
    if (this->stream_end_requested) {
	// end_stream() was called before first packet was received
	pthread_mutex_unlock(&this->state_lock);
	nupfreq = 0;
	nt_per_packet = 0;
	fpga_counts_per_sample = 0;
	fpga_count = 0;
	return false;
    }
    
    pthread_mutex_unlock(&this->state_lock);

    // After the first_packet_received flag is set, it's OK to access the fp_* fields without holding the lock.
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
    
    // Note that the arg/opaque_arg pointer is only dereferenced here, for reasons explained in a comment in make() above.
    shared_ptr<intensity_network_stream> *arg = (shared_ptr<intensity_network_stream> *) opaque_arg;
    shared_ptr<intensity_network_stream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to network_pthread_main()");

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


void intensity_network_stream::_network_thread_body()
{
    pthread_mutex_lock(&this->state_lock);

    // Advance stream state to "network_thread_started"...
    this->network_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);

    // ...and wait for it to advance to "stream_started"
    for (;;) {
	if (this->stream_end_requested) {
	    // This case can arise if end_stream() is called early
	    pthread_mutex_unlock(&this->state_lock);
	    return;
	}
	if (this->stream_started) {
	    pthread_mutex_unlock(&this->state_lock);
	    break;
	}
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    }

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

	    if (this->stream_end_requested) {
		pthread_mutex_unlock(&this->state_lock);    
                _network_flush_packets();
		return;
	    }

	    pthread_mutex_unlock(&this->state_lock);

	    // We call _add_event_counts() in a few different places in this routine, to ensure that
	    // the network thread's event counts are always regularly accumulated.
	    this->_add_event_counts(network_thread_event_subcounts);

	    cancellation_check_timestamp = curr_timestamp;
	}

	// Periodically flush packets to assembler thread (only happens if packet rate is low; normal case is that the packet_list fills first)
	if (curr_timestamp > incoming_packet_list_timestamp + constants::unassembled_ringbuf_timeout_usec) {
            _network_flush_packets();
	    incoming_packet_list_timestamp = curr_timestamp;
	}

	// Read new packet from socket (note that socket has a timeout, so this call can time out)
	uint8_t *packet_data = incoming_packet_list->data_end;

        // Read from socket, recording the sender IP & port.
        sockaddr_in sender_addr;
        int slen = sizeof(sender_addr);
        int packet_nbytes = ::recvfrom(sockfd, packet_data, constants::max_input_udp_packet_size + 1, 0,
                                       (struct sockaddr *)&sender_addr, (socklen_t *)&slen);
	// Check for error or timeout in read()
	if (packet_nbytes < 0) {
	    if ((errno == EAGAIN) || (errno == ETIMEDOUT))
		continue;  // normal timeout
	    throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
	}

        // Increment the number of packets we've received from this sender:
        uint64_t ip = ntohl(sender_addr.sin_addr.s_addr);
        uint64_t port = ntohs(sender_addr.sin_port);
        // IPv4 address in high 32 bits, port in low 16 bits
        uint64_t sender = (ip << 32) | port;
        network_thread_perhost_packets[sender]++;

	event_subcounts[event_type::byte_received] += packet_nbytes;
	event_subcounts[event_type::packet_received]++;

	// If we receive a special "short" packet (length 24), it indicates end-of-stream.
	if (_unlikely(packet_nbytes == 24)) {
	    event_subcounts[event_type::packet_end_of_stream]++;
	    if (ini_params.accept_end_of_stream_packets)
		return;   // triggers shutdown of entire stream
	    continue;
	}

	// Special logic for processing first packet, which triggers some initializations.
	if (is_first_packet) {
	    intensity_packet packet;

	    if (!packet.decode(packet_data, packet_nbytes)) {
		event_subcounts[event_type::packet_bad]++;
		continue;
	    }

	    this->fp_nupfreq = packet.nupfreq;
	    this->fp_nt_per_packet = packet.ntsamp;
	    this->fp_fpga_counts_per_sample = packet.fpga_counts_per_sample;
	    this->fp_fpga_count = packet.fpga_count;
	    is_first_packet = false;

	    // This will unblock the assembler thread.  (After being spawned, the assembler thread
	    // waits for the first_packet_received flag before doing anything else.)
	    pthread_mutex_lock(&this->state_lock);
	    this->first_packet_received = true;
	    pthread_cond_broadcast(&this->cond_state_changed);
	    pthread_mutex_unlock(&this->state_lock);
	}

	// The incoming_packet_list is timestamped with the arrival time of its first packet.
	if (incoming_packet_list->curr_npackets == 0)
	    incoming_packet_list_timestamp = curr_timestamp;

	incoming_packet_list->add_packet(packet_nbytes);

	if (incoming_packet_list->is_full)
            _network_flush_packets();
    }
}

// This gets called from the network thread to flush packets to the assembler
// threads.
void intensity_network_stream::_network_flush_packets() {
    this->_put_unassembled_packets();
    this->_add_event_counts(network_thread_event_subcounts);

    // Update the "perhost_packets" counter from "network_thread_perhost_packets"
    pthread_mutex_lock(&this->event_lock);
    for (auto it = network_thread_perhost_packets.begin();
         it != network_thread_perhost_packets.end(); it++) {
        perhost_packets[it->first] = it->second;
    }
    pthread_mutex_unlock(&this->event_lock);
}

// This gets called when the network thread exits (on all exit paths).
void intensity_network_stream::_network_thread_exit()
{
    // This just sets the stream_end_requested flag, if it hasn't been set already.
    this->end_stream();
    
    // Flush any pending packets to assembler thread.
    this->_put_unassembled_packets();
    
    // Make sure all event counts are accumulated.
    this->_add_event_counts(network_thread_event_subcounts);

    // Set end-of-stream flag in the unassembled_ringbuf, so that the assembler knows there are no more packets.
    unassembled_ringbuf->end_stream();

    // Make sure socket is closed.
    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }
}


void intensity_network_stream::_put_unassembled_packets()
{
    int npackets = incoming_packet_list->curr_npackets;

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
}


// -------------------------------------------------------------------------------------------------
//
// assembler thread


// static member function
void *intensity_network_stream::assembler_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to assembler_pthread_main()");

    // Note that the arg/opaque_arg pointer is only dereferenced here, for reasons explained in a comment in make() above.
    shared_ptr<intensity_network_stream> *arg = (shared_ptr<intensity_network_stream> *) opaque_arg;
    shared_ptr<intensity_network_stream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to assembler_thread_main()");

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


void intensity_network_stream::_assembler_thread_body()
{
    // Set the assembler_thread_started flag...
    pthread_mutex_lock(&this->state_lock);
    this->assembler_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);

    // ... and wait for first_packet_received flag to be set
    for (;;) {
   	if (this->stream_end_requested) {
	    // This case can arise if end_stream() is called early
	    pthread_mutex_unlock(&this->state_lock);
	    return;
	}
	if (this->first_packet_received) {
	    pthread_mutex_unlock(&this->state_lock);
	    break;
	}
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    }

    uint16_t nupfreq = this->fp_nupfreq;
    uint16_t nt_per_packet = this->fp_nt_per_packet;
    uint16_t fpga_counts_per_sample = fp_fpga_counts_per_sample;

    // Allocate assemblers, and set the assemblers_initialized flag.  (Note: it would be simpler to do this
    // in the network thread, since we could coalesce the 'first_packet_received' and 'assemblers_initialized'
    // states into a single state, but this led to transient packet loss which was problematic for unit tests.)

    this->assemblers.resize(nassemblers);
    for (int ix = 0; ix < nassemblers; ix++)
	assemblers[ix] = make_unique<assembled_chunk_ringbuf>(ini_params, ini_params.beam_ids[ix], nupfreq, nt_per_packet, fpga_counts_per_sample, this->fp_fpga_count);

    pthread_mutex_lock(&this->state_lock);
    this->assemblers_initialized = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);
    
    auto packet_list = make_unique<udp_packet_list> (constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list);
    int64_t *event_subcounts = &this->assembler_thread_event_subcounts[0];

    // Main packet loop

    while (unassembled_ringbuf->get_packet_list(packet_list)) {
	for (int ipacket = 0; ipacket < packet_list->curr_npackets; ipacket++) {
            uint8_t *packet_data = packet_list->get_packet_data(ipacket);
            int packet_nbytes = packet_list->get_packet_nbytes(ipacket);
	    intensity_packet packet;

	    if (!packet.decode(packet_data, packet_nbytes)) {
		event_subcounts[event_type::packet_bad]++;
		continue;
	    }

	    bool mismatch = (packet.nupfreq != nupfreq) || (packet.ntsamp != nt_per_packet) || (packet.fpga_counts_per_sample != fpga_counts_per_sample);

	    if (_unlikely(mismatch)) {
		event_subcounts[event_type::first_packet_mismatch]++;
		if (ini_params.throw_exception_on_first_packet_mismatch)
		    throw runtime_error("ch_frb_io: first-packet mismatch occurred and stream was constructed with 'throw_exception_on_first_packet_mismatch' flag");
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
	    // only beam index zero.  This scheme avoids the overhead of copying the packet.
	    
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
			if (ini_params.throw_exception_on_beam_id_mismatch)
			    throw runtime_error("ch_frb_io: beam_id mismatch occurred and stream was constructed with 'throw_exception_on_beam_id_mismatch' flag");
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

	// We accumulate event counts once per udp_packet_list.
	this->_add_event_counts(assembler_thread_event_subcounts);
    }
}


// Called whenever the assembler thread exits (on all exit paths)
void intensity_network_stream::_assembler_thread_exit()
{
    this->end_stream();
    unassembled_ringbuf->end_stream();

    // Use assemblers.size() instead of 'nassemblers', to correctly handle the case where end_stream()
    // gets called early, and the assemblers never get allocated.  (We also test for empty pointers, to
    // handle a corner case where an exception is thrown in the middle of the ring buffer allocation.)

    for (unsigned int i = 0; i < assemblers.size(); i++) {
	if (assemblers[i])
	    assemblers[i]->end_stream(&assembler_thread_event_subcounts[0]);
    }

    // Make sure all event counts are accumulated.
    this->_add_event_counts(assembler_thread_event_subcounts);

    // The rest of this routine just prints a summary.

    vector<int64_t> counts = this->get_event_counts();

    stringstream ss;
    ss << "ch_frb_io: assembler thread exiting\n"
       << "    bytes received (GB): " << (1.0e-9 * counts[event_type::byte_received]) << "\n"
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
