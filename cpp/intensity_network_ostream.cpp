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
// class intensity_network_ostream



// static member function (de facto constructor)
auto intensity_network_ostream::make(const std::string &dstname_, const std::vector<int> &beam_ids_, 
				     const std::vector<int> &coarse_freq_ids_, int nupfreq_, int nt_per_chunk_,
				     int nfreq_coarse_per_packet_, int nt_per_packet_, int fpga_counts_per_sample_,
				     float wt_cutoff_, double target_gbps_) 
    -> std::shared_ptr<intensity_network_ostream>
{
    auto p = new intensity_network_ostream(dstname_, beam_ids_, coarse_freq_ids_, nupfreq_, 
					   nt_per_chunk_, nfreq_coarse_per_packet_, nt_per_packet_, 
					   fpga_counts_per_sample_, wt_cutoff_, target_gbps_);
    
    shared_ptr<intensity_network_ostream> ret(p);

    ret->_open_socket();

    int err = pthread_create(&ret->network_thread, NULL, intensity_network_ostream::network_pthread_main, (void *) &ret);
    if (err < 0)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_ostream constructor: ") + strerror(errno));

    pthread_mutex_lock(&ret->state_lock);
    while (!ret->network_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->state_lock);
    pthread_mutex_unlock(&ret->state_lock);

    return ret;
}

    
intensity_network_ostream::intensity_network_ostream(const std::string &dstname_, const std::vector<int> &beam_ids_, 
						     const std::vector<int> &coarse_freq_ids_, int nupfreq_, int nt_per_chunk_,
						     int nfreq_coarse_per_packet_, int nt_per_packet_, int fpga_counts_per_sample_,
						     float wt_cutoff_, double target_gbps_) :
    dstname(dstname_),
    nbeams(beam_ids_.size()),
    nfreq_coarse_per_chunk(coarse_freq_ids_.size()),
    nfreq_coarse_per_packet(nfreq_coarse_per_packet_),
    nupfreq(nupfreq_),
    nt_per_chunk(nt_per_chunk_),
    nt_per_packet(nt_per_packet_),
    fpga_counts_per_sample(fpga_counts_per_sample_),
    nbytes_per_packet(packet_size(nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet)),
    npackets_per_chunk((nfreq_coarse_per_chunk / nfreq_coarse_per_packet) * (nt_per_chunk / nt_per_packet)),
    nbytes_per_chunk(nbytes_per_packet * npackets_per_chunk),
    wt_cutoff(wt_cutoff_),
    target_gbps(target_gbps_),
    beam_ids(beam_ids_),
    coarse_freq_ids(coarse_freq_ids_)
{
    // Tons of argument checking.

    if ((beam_ids.size() == 0) || (beam_ids.size() >= 65536))
	throw runtime_error("chime intensity_network_ostream constructor: beam_ids vector is empty or too large");

    if ((coarse_freq_ids.size() == 0) || (coarse_freq_ids.size() >= 65536))
	throw runtime_error("chime intensity_network_ostream constructor: coarse_freq_ids vector is empty or too large");
    if (nfreq_coarse_per_packet <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nfreq_per_packet > 0");
    if (nfreq_coarse_per_chunk % nfreq_coarse_per_packet != 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nfreq_per_chunk to be a multiple of nfreq_per_packet");

    if (nt_per_chunk <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_chunk > 0");
    if (nt_per_packet <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_packet > 0");
    if (!is_power_of_two(nt_per_packet))
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_packet to be a power of two");
    if (nt_per_chunk % nt_per_packet != 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_chunk to be a multiple of nt_per_packet");

    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of nupfreq");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of fpga_counts_per_sample");
    if (wt_cutoff < 0.0)
	throw runtime_error("chime intensity_network_ostream constructor: expected wt_cutoff to be >= 0.0");
    if (target_gbps < 0.0)
	throw runtime_error("chime intensity_network_ostream constructor: expected target_gbps to be >= 0.0");

    if (nbytes_per_packet > constants::max_output_udp_packet_size)
	throw runtime_error("chime intensity_network_ostream constructor: packet size is too large, you need to decrease nfreq_per_packet or nt_per_packet");

    for (unsigned int i = 0; i < beam_ids.size(); i++) {
	if ((beam_ids[i] < 0) || (beam_ids[i] > constants::max_allowed_beam_id))
	    throw runtime_error("intensity_network_ostream constructor: bad beam_id");
	for (unsigned int j = 0; j < i; j++)
	    if (beam_ids[i] == beam_ids[j])
		throw runtime_error("intensity_network_ostream constructor: duplicate beam_id");
    }

    for (unsigned int i = 0; i < coarse_freq_ids.size(); i++) {
	if ((coarse_freq_ids[i] < 0) || (coarse_freq_ids[i] >= constants::nfreq_coarse))
	    throw runtime_error("intensity_network_ostream constructor: bad coarse_freq_id");
	for (unsigned int j = 0; j < i; j++)
	    if (coarse_freq_ids[i] == coarse_freq_ids[j])
		throw runtime_error("intensity_network_ostream constructor: duplicate coarse_freq_id");
    }

    // Parse dstname: expect string of the form HOSTNAME[:PORT]

    this->hostname = dstname;
    this->udp_port = constants::default_udp_port;

    size_t i = dstname.find(":");

    if (i != std::string::npos) {
	string portstr = dstname.substr(i+1);
	try {
	    this->udp_port = lexical_cast<uint16_t> (portstr);	    
	} catch (...) {
	    throw runtime_error("ch_frb_io: couldn't convert string '" + portstr + "' to 16-bit udp port number");
	}

	this->hostname = dstname.substr(0,i);
    }

    // Remaining initializations (except socket, which is initialized in intensity_network_ostream::_open_socket())

    this->beam_ids_16bit.resize(nbeams, 0);
    for (int i = 0; i < nbeams; i++)
	beam_ids_16bit[i] = uint16_t(beam_ids[i]);

    this->coarse_freq_ids_16bit.resize(nfreq_coarse_per_chunk, 0);
    for (int i = 0; i < nfreq_coarse_per_chunk; i++)
	coarse_freq_ids_16bit[i] = uint16_t(coarse_freq_ids[i]);
    
    xpthread_mutex_init(&this->state_lock);
    xpthread_cond_init(&this->cond_state_changed);

    int capacity = constants::output_ringbuf_capacity;
    this->ringbuf = make_unique<udp_packet_ringbuf> (capacity, npackets_per_chunk, nbytes_per_chunk);
    this->tmp_packet_list = udp_packet_list(npackets_per_chunk, nbytes_per_chunk);
}


intensity_network_ostream::~intensity_network_ostream()
{
    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }

    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&state_lock);
}


// Socket initialization factored to its own routine, rather than putting it in the constructor,
// so that the socket will always be closed if an exception is thrown somewhere.
void intensity_network_ostream::_open_socket()
{
    if (sockfd >= 0)
	throw runtime_error("double call to intensity_network_ostream::_open_socket()");

    struct sockaddr_in saddr;
    memset(&saddr, 0, sizeof(saddr));
    saddr.sin_family = AF_INET;
    saddr.sin_port = htons(this->udp_port);
    
    // FIXME need getaddrinfo() here
    int err = inet_pton(AF_INET, this->hostname.c_str(), &saddr.sin_addr);
    if (err == 0)
	throw runtime_error("ch_frb_io: couldn't resolve hostname '" + hostname + "' to an IP address: general parse error");
    if (err < 0)
	throw runtime_error("ch_frb_io: couldn't resolve hostname '" + hostname + "' to an IP address: " + strerror(errno) + "general parse error");

    this->sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0)
	throw runtime_error(string("ch_frb_io: couldn't create udp socket: ") + strerror(errno));

    int socket_bufsize = constants::send_socket_bufsize;
    err = setsockopt(sockfd, SOL_SOCKET, SO_SNDBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_SNDBUF) failed: ") + strerror(errno));
    
    // Note: bind() not called, so source port number of outgoing packets will be arbitrarily assigned

    if (connect(sockfd, reinterpret_cast<struct sockaddr *> (&saddr), sizeof(saddr)) < 0)
	throw runtime_error("ch_frb_io: couldn't connect udp socket to dstname '" + dstname + "': " + strerror(errno));
}


// The 'intensity' and 'weights' arrays have shapes (nbeams, nfreq_coarse_per_chunk, nupfreq, nt_per_chunk)
void intensity_network_ostream::send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count, bool is_blocking)
{
    if (fpga_count % (fpga_counts_per_sample * nt_per_packet) != 0)
	throw runtime_error("intensity_network_ostream::send_chunk(): fpga count must be divisible by (fpga_counts_per_sample * nt_per_packet)");
    if (tmp_packet_list.curr_npackets > 0)
	throw runtime_error("intensity_network_ostream::send_chunk(): internal error: tmp_packet_list nonempty?!");

    intensity_packet packet;
    
    // Some intensity_packet fields are packet-independent; these are initialized here.
    packet.protocol_version = 1;
    packet.data_nbytes = nbeams * nfreq_coarse_per_packet * nupfreq * nt_per_packet;
    packet.fpga_counts_per_sample = fpga_counts_per_sample;
    packet.nbeams = nbeams;
    packet.nfreq_coarse = nfreq_coarse_per_packet;
    packet.nupfreq = nupfreq;
    packet.ntsamp = nt_per_packet;
    packet.beam_ids = &beam_ids_16bit[0];

    // The number of packets per chunk is (nf_outer * nt_outer)
    int beam_stride = nfreq_coarse_per_chunk * nupfreq * stride;
    int nf_outer = nfreq_coarse_per_chunk / nfreq_coarse_per_packet;
    int nt_outer = nt_per_chunk / nt_per_packet;

    // Outer loop over packets in chunk
    for (int if_outer = 0; if_outer < nf_outer; if_outer++) {
	for (int it_outer = 0; it_outer < nt_outer; it_outer++) {
	    int data_offset = (if_outer * nfreq_coarse_per_packet * nupfreq * stride) + (it_outer * nt_per_packet);

	    // Some intensity_packet fields are packet-dependent; these are initialized here.
	    packet.freq_ids = &coarse_freq_ids_16bit[if_outer * nfreq_coarse_per_packet];
	    packet.fpga_count = fpga_count + it_outer * nt_per_packet * fpga_counts_per_sample;

	    packet.encode(tmp_packet_list.data_end, intensity + data_offset, weights + data_offset, beam_stride, stride, wt_cutoff);
	    tmp_packet_list.add_packet(nbytes_per_packet);
	}
    }

    // FIXME should count dropped packets here
    // FIXME should have boolean drops_allowed
    ringbuf->put_packet_list(tmp_packet_list, is_blocking);
}


void intensity_network_ostream::end_stream(bool join_network_thread)
{
    ringbuf->end_stream();

    if (!join_network_thread)
	return;

    pthread_mutex_lock(&this->state_lock);

    if (network_thread_joined) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: attempt to join ostream output thread twice");
    }
    
    pthread_mutex_unlock(&this->state_lock);
    pthread_join(network_thread, NULL);
}


// -------------------------------------------------------------------------------------------------
//
// Network write thread


// static member function
void *intensity_network_ostream::network_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_pthread_main()");

    shared_ptr<intensity_network_ostream> *arg = (shared_ptr<intensity_network_ostream> *) opaque_arg;
    shared_ptr<intensity_network_ostream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to network_pthread_main()");

    stream->_network_thread_start();

    try {
	stream->_network_thread_body();
    } catch (...) {
	stream->end_stream(false);   // "false" means "don't join threads" (would deadlock otherwise!)
	throw;
    }

    stream->end_stream(false);   // "false" has same meaning as above
    return NULL;
}


void intensity_network_ostream::_network_thread_start()
{
    pthread_mutex_lock(&state_lock);
    network_thread_started = true;
    pthread_cond_broadcast(&cond_state_changed);
    pthread_mutex_unlock(&state_lock);
}


void intensity_network_ostream::_network_thread_body()
{
    udp_packet_list packet_list(npackets_per_chunk, nbytes_per_chunk);
    int last_packet_nbytes = 0;

    // to be initialized when first packet is sent
    struct timeval tv_ini;
    
    // Loop over packet_lists
    for (;;) {
	if (!ringbuf->get_packet_list(packet_list))
	    break;   // end of stream reached (probably normal termination)
	
	// Loop over packets
	for (int ipacket = 0; ipacket < packet_list.curr_npackets; ipacket++) {
	    const uint8_t *packet = packet_list.data_start + packet_list.packet_offsets[ipacket];
	    const int packet_nbytes = packet_list.packet_offsets[ipacket+1] - packet_list.packet_offsets[ipacket];

	    if (npackets_sent == 0)
		tv_ini = xgettimeofday();

	    int64_t last_timestamp = this->curr_timestamp;
	    curr_timestamp = usec_between(tv_ini, xgettimeofday());
	    
	    // Throttling logic: compare actual bandwidth to 'target_gbps' and sleep if necessary.
	    if ((target_gbps > 0.0) && (npackets_sent > 0)) {
		int64_t target_timestamp = last_timestamp + int64_t(8.0e-3 * last_packet_nbytes / target_gbps);		
		if (curr_timestamp < target_timestamp) {
		    xusleep(target_timestamp - curr_timestamp);
		    curr_timestamp = target_timestamp;
		}
	    }

	    ssize_t n = send(this->sockfd, packet, packet_nbytes, 0);
	    
	    if (n < 0)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() failed: ") + strerror(errno));
	    if (n != packet_nbytes)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() sent ") + to_string(n) + "/" + to_string(packet_nbytes) + " bytes?!");

	    last_packet_nbytes = packet_nbytes;
	    this->nbytes_sent += packet_nbytes;
	    this->npackets_sent++;
	}
    }
    
    this->_announce_end_of_stream();
}


void intensity_network_ostream::_announce_end_of_stream()
{
    // FIXME temporary hack.  For testing, it is convenient to have a way of ending the stream. 
    // We make the rule that a packet with nbeams = nfreq_coarse = nupfreq = ntsamp = 0 means "end of stream".
    // Later this will be replaced by something better!  
    //
    // Since UDP doesn't guarantee delivery, we have no way to ensure that the end-of-stream packet 
    // reaches the other side, but we'll make a best effort by sending 10 packets separated by 0.1 sec.

    cerr << "ch_frb_io: network output thread sending end-of-stream packets\n";

    for (int ipacket = 0; ipacket < 10; ipacket++) {
	vector<uint8_t> packet(24, uint8_t(0));
	*((uint32_t *) &packet[0]) = uint32_t(1);  // protocol number

	ssize_t n = send(this->sockfd, &packet[0], packet.size(), 0);

	if (n == (ssize_t)packet.size()) {
	    usleep(100000);  // 10^5 microseconds
	    continue;
	}

	// Emit warning if we fail on the first packet.  No warning emitted if subsequent packets
	// fail, since a likely explanation is that we're running over the loopback interface and
	// the receiving socket has been closed.

	if (ipacket == 0)
	    cerr << "warning: end-of-stream packets failed to send\n";

	break;
    }

    // End-of-stream summary info

    stringstream ss;
    ss << "ch_frb_io: network output thread exiting, npackets_sent=" << npackets_sent;
    if (npackets_sent >= 2)
	ss << ", gbps=" << (8.0e-3 * nbytes_sent / double(curr_timestamp));
    if (target_gbps > 0.0)
	ss << ", target_gbps=" << target_gbps;
    ss << "\n";

    cerr << ss.str();
}


}  // namespace ch_frb_io
