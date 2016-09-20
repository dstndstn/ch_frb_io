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

using event_counts = intensity_network_stream::event_counts;

void event_counts::clear()
{
    this->num_bad_packets = 0;
    this->num_good_packets = 0;
    this->num_beam_id_mismatches = 0;
    this->num_first_packet_mismatches = 0;
}

event_counts &event_counts::operator+=(const event_counts &x)
{
    this->num_bad_packets += x.num_bad_packets;
    this->num_good_packets += x.num_good_packets;
    this->num_beam_id_mismatches += x.num_beam_id_mismatches;
    this->num_first_packet_mismatches += x.num_first_packet_mismatches;
    return *this;
}


// -------------------------------------------------------------------------------------------------
//
// intensity_network_stream


// Static member function (de facto constructor)
shared_ptr<intensity_network_stream> intensity_network_stream::make(const vector<shared_ptr<intensity_beam_assembler> > &assemblers, int udp_port)
{
    intensity_network_stream *retp = new intensity_network_stream(assemblers, udp_port);
    shared_ptr<intensity_network_stream> ret(retp);

    ret->_open_socket();

    // Spawn assembler thread
    int err = pthread_create(&ret->assembler_thread, NULL, assembler_pthread_main, (void *) &ret);
    if (err)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_stream constructor: ") + strerror(errno));
    
    // Wait for assembler thread to start
    pthread_mutex_lock(&ret->lock);
    while (!ret->assembler_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->lock);
    pthread_mutex_unlock(&ret->lock);    

    // Spawn network thread
    err = pthread_create(&ret->network_thread, NULL, network_pthread_main, (void *) &ret);
    if (err)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_stream constructor: ") + strerror(errno));
    
    // Wait for network thread to start
    pthread_mutex_lock(&ret->lock);
    while (!ret->network_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->lock);
    pthread_mutex_unlock(&ret->lock);    

    return ret;
}


intensity_network_stream::intensity_network_stream(const vector<shared_ptr<intensity_beam_assembler> > &assemblers_, int udp_port_) :
    udp_port(udp_port_),
    incoming_packet_list(constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list),
    assemblers(assemblers_),
    nassemblers(assemblers_.size())
{
    // Argument checking

    if (nassemblers == 0)
	throw runtime_error("ch_frb_io: empty assembler list passed to intensity_network_stream constructor");

    for (int i = 0; i < nassemblers; i++) {
	if (!assemblers_[i])
	    throw runtime_error("ch_frb_io: empty assembler pointer passed to intensity_network_stream constructor");

	for (int j = 0; j < i; j++)
	    if (assemblers_[i]->beam_id == assemblers_[j]->beam_id)
		throw runtime_error("ch_frb_io: assembler list passed to intensity_network_stream constructor contains duplicate beam_ids");
    }

    if ((udp_port <= 0) || (udp_port >= 65536))
	throw runtime_error("ch_frb_io: intensity_network_stream constructor: bad udp port " + to_string(udp_port));

    // All initializations except the socket (which is initialized in _open_socket())

    int capacity = constants::unassembled_ringbuf_capacity;
    int max_npackets = constants::max_unassembled_packets_per_list;
    int max_nbytes = constants::max_unassembled_nbytes_per_list;
    string dropstr = "warning: assembler thread running too slow, dropping packets";
    bool drops_allowed = true;    // FIXME
    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (capacity, max_npackets, max_nbytes, dropstr, drops_allowed);

    this->assembler_beam_ids.resize(nassemblers, -1);
    for (int i = 0; i < nassemblers; i++)
	assembler_beam_ids[i] = assemblers[i]->beam_id;

    pthread_mutex_init(&lock, NULL);
    pthread_cond_init(&cond_state_changed, NULL);
}


intensity_network_stream::~intensity_network_stream()
{
    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&lock);

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


// Just sets the stream_ended flag and returns.
// The network thread will see that this flag has been set, flush packets to the assembler thread, call unassembled_ringbuf->end_stream(), and exit.
// The assembler thread will see that the ringbuf has ended, flush buffers to the processing threads, and exit.
void intensity_network_stream::end_stream()
{
    pthread_mutex_lock(&this->lock);
    this->stream_started = true;
    this->stream_ended = true;    
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->lock);
}


void intensity_network_stream::join_threads()
{
    pthread_mutex_lock(&this->lock);
    
    if (!stream_started) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::join_network_thread() was called with no prior call to start_stream()");
    }

    if (join_called) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: double call to intensity_network_stream::join_network_thread()");
    }

    this->join_called = true;
    pthread_mutex_unlock(&this->lock);

    pthread_join(network_thread, NULL);
    pthread_join(assembler_thread, NULL);
}


intensity_network_stream::event_counts intensity_network_stream::get_event_counts() const
{
    pthread_mutex_lock(&this->lock);
    event_counts ret = this->curr_counts;
    pthread_mutex_unlock(&this->lock);    
    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


// FIXME improve?
inline bool is_end_of_stream_packet(const uint8_t *packet, int packet_nbytes)
{
    return packet_nbytes == 24;
}


// static member function
void *intensity_network_stream::network_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_pthread_main()");

    shared_ptr<intensity_network_stream> *arg = (shared_ptr<intensity_network_stream> *) opaque_arg;
    shared_ptr<intensity_network_stream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to network_thread_main()");

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

    pthread_mutex_lock(&this->lock);

    this->network_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);

    for (;;) {
	if (this->stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return;
	}
	if (this->stream_started) {
	    pthread_mutex_unlock(&this->lock);
	    break;
	}
	pthread_cond_wait(&this->cond_state_changed, &this->lock);
    }
}


void intensity_network_stream::_network_thread_body()
{
    // Start listening on socket 

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));
	
    server_address.sin_family = AF_INET;
    inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
    server_address.sin_port = htons(udp_port);

    int err = ::bind(sockfd, (struct sockaddr *) &server_address, sizeof(server_address));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: bind() failed: ") + strerror(errno));

    cerr << ("ch_frb_io: listening for packets on port " + to_string(udp_port) + "\n");

    // Main packet loop

    struct timeval tv_ini = xgettimeofday();
    uint64_t packet_list_start_timestamp = 0;
    uint64_t cancellation_check_timestamp = 0;

    for (;;) {
	// Read new packet from socket (note that socket has a timeout, so this call can time out)
	uint8_t *packet = incoming_packet_list.data_end;
	int packet_nbytes = read(sockfd, packet, constants::max_input_udp_packet_size + 1);

	// All timestamps are in microseconds relative to tv_ini.
	uint64_t curr_timestamp = usec_between(tv_ini, xgettimeofday());
	
	// Periodically check whether stream has been cancelled by end_stream().
	if (curr_timestamp > cancellation_check_timestamp + constants::stream_cancellation_latency_usec) {
	    pthread_mutex_lock(&this->lock);

	    if (!this->stream_ended) {
		pthread_mutex_unlock(&this->lock);    
		return;
	    }

	    pthread_mutex_unlock(&this->lock);
	    cancellation_check_timestamp = curr_timestamp;
	}

	// Periodically flush packets to assembler thread (only happens if packet rate is low; normal case is that the packet_list fills first)
	if (curr_timestamp > packet_list_start_timestamp + constants::unassembled_ringbuf_timeout_usec)
	    this->_put_unassembled_packets();

	// Check for error or timeout in read()
	if (packet_nbytes < 0) {
	    if ((errno == EAGAIN) || (errno == ETIMEDOUT))
		continue;  // normal timeout
	    throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
	}

	// If we receive a special "short" packet (length 24), it indicates end-of-stream.
	// FIXME is this a temporary kludge or something which should be documented in the packet protocol?
	if (is_end_of_stream_packet(packet, packet_nbytes))
	    return;

	// If we get here, then the packet will be sent to the assembler thread

	if (incoming_packet_list.curr_npackets == 0)
	    packet_list_start_timestamp = curr_timestamp;

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

    this->_put_unassembled_packets();
    unassembled_ringbuf->end_stream();
}


void intensity_network_stream::_put_unassembled_packets()
{
    if (incoming_packet_list.curr_npackets == 0)
	return;
    
    if (!unassembled_ringbuf->producer_put_packet_list(incoming_packet_list, false))
	throw runtime_error("ch_frb_io: internal error: couldn't write packets to unassembled_ringbuf");

    // FIXME maintain counts of packets_received, packets_dropped
    // FIXME throw exception if drops_allowed=false and packets were dropped
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
    pthread_mutex_lock(&this->lock);
    this->assembler_thread_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->lock);
}


void intensity_network_stream::_assembler_thread_body()
{
    udp_packet_list packet_list(constants::max_unassembled_packets_per_list, constants::max_unassembled_nbytes_per_list);

    while (unassembled_ringbuf->consumer_get_packet_list(packet_list)) {
	this->_tmp_counts.clear();

	for (int ipacket = 0; ipacket < packet_list.curr_npackets; ipacket++) {
            uint8_t *packet_data = packet_list.data_start + packet_list.packet_offsets[ipacket];
            int packet_nbytes = packet_list.packet_offsets[ipacket+1] - packet_list.packet_offsets[ipacket];
	    
	    this->_assemble_packet(packet_data, packet_nbytes);
	}

	pthread_mutex_lock(&this->lock);
	this->curr_counts += _tmp_counts;
	pthread_mutex_unlock(&this->lock);
    }
}


void intensity_network_stream::_assembler_thread_exit()
{
    event_counts counts = this->get_event_counts();

    stringstream ss;
    ss << "ch_frb_io: network input thread exiting:"
       << " good_packets=" << counts.num_good_packets 
       << ", bad_packets=" << counts.num_bad_packets 
       << ", beam_id_mismatches=" << counts.num_beam_id_mismatches 
       << "\n";

    cerr << ss.str().c_str();

    unassembled_ringbuf->end_stream();

    for (int i = 0; i < nassemblers; i++)
	assemblers[i]->_end_stream();
}


void intensity_network_stream::_assemble_packet(const uint8_t *packet_data, int packet_nbytes)
{
    intensity_packet packet;

    if (!packet.read(packet_data, packet_nbytes)) {
	this->_tmp_counts.num_bad_packets++;
	return;
    }

    // The following checks are not included in intensity_packet::read().

    if (!is_power_of_two(packet.ntsamp)) {
	this->_tmp_counts.num_bad_packets++;
	return;
    }

    for (int i = 0; i < packet.nfreq_coarse; i++) {
	if (_unlikely(packet.freq_ids[i] >= constants::nfreq_coarse)) {
	    this->_tmp_counts.num_bad_packets++;
	    return;
	}
    }

    if (this->expected_params_initialized) {
	bool mismatch = ((packet.nupfreq != expected_nupfreq) || 
			 (packet.ntsamp != expected_nt_per_packet) ||
			 (packet.fpga_counts_per_sample != expected_fpga_counts_per_sample));
	
	if (_unlikely(mismatch)) {
	    this->_tmp_counts.num_first_packet_mismatches++;
	    return;
	}
    }
    else {
	if (_unlikely(packet.fpga_counts_per_sample == 0)) {
	    this->_tmp_counts.num_bad_packets++;
	    return;
	}

	this->expected_nupfreq = packet.nupfreq;
	this->expected_nt_per_packet = packet.ntsamp;
	this->expected_fpga_counts_per_sample = packet.fpga_counts_per_sample;
	this->expected_params_initialized = true;
    }

    // Note conversions to uint64_t, to prevent integer overflow
    uint64_t fpga_counts_per_packet = uint64_t(packet.fpga_counts_per_sample) * uint64_t(packet.ntsamp);
    if (_unlikely(packet.fpga_count % fpga_counts_per_packet != 0)) {
	this->_tmp_counts.num_bad_packets++;
	return;
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

    this->_tmp_counts.num_good_packets++;

    int nbeams = packet.nbeams;
    int nfreq_coarse = packet.nfreq_coarse;
    int new_data_nbytes = nfreq_coarse * packet.nupfreq * packet.ntsamp;
    int *assembler_ids = &assembler_beam_ids[0];  // bare pointer for speed

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
	    if (assembler_ids[assembler_ix] == packet_id) {
		// Match found
		assemblers[assembler_ix]->_put_unassembled_packet(packet);
		break;
	    }

	    assembler_ix++;

	    if (assembler_ix >= nassemblers) {
		// No match found
		this->_tmp_counts.num_beam_id_mismatches++;
		break;
	    }
	}

	// Danger zone: we do some pointer arithmetic, to modify the packet so that it now
	// corresponds to a new subset of the original packet, corresponding to beam index (ibeam+1).

	packet.beam_ids += 1;
	packet.scales += nfreq_coarse;
	packet.offsets += nfreq_coarse;
	packet.data += new_data_nbytes;
    }
}


}  // namespace ch_frb_io
