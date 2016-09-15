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
static void network_thread_main2(intensity_network_ostream *stream);



// -------------------------------------------------------------------------------------------------
//
// Helper functions for encoding L0_L1 packets.


//
// Encodes one "row" of the L0_L1 packet.  A row is a subset of the data with a common (offset, scale).
// In the current packet protocol this corresponds to all (upfreq, time) pairs for a fixed (beam, freq).
//
// The 'intensity' array is a contiguous array of length 'rowlen'.
// The 'mask' array is a contiguous array of length 'rowlen', which is assumed to contain zeros or ones.
//
// We currently assign the scale by computing the unclipped variance and saturating at 5 sigma.
// FIXME this could be improved by clipping 5 sigma outliers and iteratively computing the variance.
// 
inline void encode_packet_row(int rowlen, float *scalep, float *offsetp, uint8_t *datap, const float *intensity, const float *mask)
{
    double acc0 = 0.0;
    double acc1 = 0.0;
    double acc2 = 0.0;

    for (int i = 0; i < rowlen; i++) {
	acc0 += mask[i];
	acc1 += mask[i] * intensity[i];
	acc2 += mask[i] * intensity[i] * intensity[i];
    }

    if (acc0 <= 0.0) {
	*scalep = 1.0;
	*offsetp = 0.0;
	memset(datap, 0, rowlen);
	return;
    }

    double mean = acc1/acc0;
    double var = acc2/acc0 - mean*mean;
    
    var = max(var, double(1.0e-10*mean*mean));

    double scale = sqrt(var) / 25.;
    double offset = -128.*scale + mean;   // 0x80 -> mean

    *scalep = scale;
    *offsetp = offset;

    for (int i = 0; i < rowlen; i++) {
	float t = (intensity[i] - offset) / scale;
	t *= mask[i];           // so that masked samples get encoded as 0x00
	t = min(t, float(255.));
	t = max(t, float(0.));
	datap[i] = int(t+0.5);  // round to nearest integer
    }
}

//
// Note: encode_packet() doesn't do much argument checking.
//
// Some checks that the caller is responsible for:
//   - nbeam, nfreq, nupfreq, ntsamp should all be >=1 and their product should be <= 8900 (udp mtu)
//   - fpga_counts_per_sample should fit into 16 bits
//   - fpga_count should be divisible by fpga_counts_per_sample
//
inline void encode_packet(int nbeam, int nfreq, int nupfreq, int ntsamp, 
			  uint16_t fpga_counts_per_sample, uint64_t fpga_count,
			  uint8_t *out, const uint16_t *ibeam, const uint16_t *ifreq, 
			  const float *intensity, const float *mask)
{
    int data_nbytes = nbeam * nfreq * nupfreq * ntsamp;

    // Write 24-byte header

    *((uint32_t *) out) = 1;                            //  uint32_t  protocol
    *((int16_t *) (out+4)) = data_nbytes;               //  int16_t   data_nbytes
    *((uint16_t *) (out+6)) = fpga_counts_per_sample;   //  uint16_t  fpga_counts_per_sample
    *((uint64_t *) (out+8)) = fpga_count;               //  uint64_t  fpga_count
    *((uint16_t *) (out+16)) = uint16_t(nbeam);         //  uint16_t  nbeam
    *((uint16_t *) (out+18)) = uint16_t(nfreq);         //  uint16_t  nfreq
    *((uint16_t *) (out+20)) = uint16_t(nupfreq);       //  uint16_t  nupfreq
    *((uint16_t *) (out+22)) = uint16_t(ntsamp);        //  uint16_t  ntsamp

    // Write ibeam, ifreq arrays
    // Byte count = (2*nbeam + 2*nfreq)

    memcpy(out + 24, ibeam, 2*nbeam);
    memcpy(out + 24 + 2*nbeam, ifreq, 2*nfreq);

    // Write offset, scale, data
    // Byte count = (8*nbeam*nfreq) + (nbeam*nfreq*nupfreq*nt)

    static_assert(sizeof(float)==4, "the logic below assumes sizeof(float)==4");

    // In this routine, a "row" is a (beam, freq_coarse) pair.
    int nrows = nbeam * nfreq;
    int rowlen = nupfreq * ntsamp;

    float *scale0 = (float *) (out + 24 + 2*nbeam + 2*nfreq);
    float *offset0 = (float *) (out + 24 + 2*nbeam + 2*nfreq + 4*nrows);
    uint8_t *data0 = out + 24 + 2*nbeam + 2*nfreq + 8*nrows;

    for (int irow = 0; irow < nrows; irow++)
	encode_packet_row(rowlen, scale0+irow, offset0+irow, data0 + irow*rowlen, intensity + irow*rowlen, mask + irow*rowlen);
}


// -------------------------------------------------------------------------------------------------
//
// Helper functions for intensity_network_ostream constructor


static int make_socket_from_dstname(string dstname)
{
    // Parse dstname: expect string of the form HOSTNAME[:PORT]

    size_t i = dstname.find(":");
    uint16_t port = constants::default_udp_port;

    if (i != std::string::npos) {
	string portstr = dstname.substr(i+1);
	try {
	    port = lexical_cast<uint16_t> (portstr);	    
	} catch (...) {
	    throw runtime_error("ch_frb_io: couldn't convert string '" + portstr + "' to 16-bit udp port number");
	}

	// remove port number
	dstname = dstname.substr(0,i);
    }

    struct sockaddr_in saddr;
    memset(&saddr, 0, sizeof(saddr));
    saddr.sin_family = AF_INET;
    saddr.sin_port = htons(port);
    
    // FIXME need getaddrinfo() here
    int err = inet_pton(AF_INET, dstname.c_str(), &saddr.sin_addr);
    if (err == 0)
	throw runtime_error("ch_frb_io: couldn't resolve hostname '" + dstname + "' to an IP address: general parse error");
    if (err < 0)
	throw runtime_error("ch_frb_io: couldn't resolve hostname '" + dstname + "' to an IP address: " + strerror(errno) + "general parse error");

    int sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0)
	throw runtime_error(string("ch_frb_io: couldn't create udp socket: ") + strerror(errno));

    int socket_bufsize = constants::send_socket_bufsize;
    err = setsockopt(sockfd, SOL_SOCKET, SO_SNDBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_SNDBUF) failed: ") + strerror(errno));
    
    // Note: bind() not called, so source port number of outgoing packets will be arbitrarily assigned

    if (connect(sockfd, reinterpret_cast<struct sockaddr *> (&saddr), sizeof(saddr)) < 0) {
	close(sockfd);
	throw runtime_error("ch_frb_io: couldn't connect udp socket to dstname '" + dstname + "': " + strerror(errno));
    }

    return sockfd;
}


inline vector<uint16_t> _init_ivec(const vector<int> &v, const string &objname)
{
    if (v.size() == 0)
	throw runtime_error("intensity_network_ostream constructor: " + objname + " index array is empty");	
    if (v.size() >= 65536)
	throw runtime_error("intensity_network_ostream constructor: " + objname + " index array length is >= 2^16?!");

    int nelts = v.size();
    vector<uint16_t> ret(nelts);

    for (int i = 0; i < nelts; i++) {
	if ((v[i] < 0) || (v[i] >= 65536))
	    throw runtime_error("intensity_network_ostream constructor: " + objname + " indices must be between 0 and 2^16");
	ret[i] = uint16_t(v[i]);
    }

    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// intensity_network_ostream

    
intensity_network_ostream::intensity_network_ostream(const std::string &dstname, const std::vector<int> &ibeam_, 
						     const std::vector<int> &ifreq_chunk_, int nupfreq_, int nt_per_chunk_,
						     int nfreq_per_packet_, int nt_per_packet_, int fpga_counts_per_sample_,
						     float wt_cutoff_, double target_gbps_) :
    sockfd(make_socket_from_dstname(dstname)),
    nbeam(ibeam_.size()),
    nupfreq(nupfreq_),
    nfreq_per_chunk(ifreq_chunk_.size()),
    nfreq_per_packet(nfreq_per_packet_),
    nt_per_chunk(nt_per_chunk_),
    nt_per_packet(nt_per_packet_),
    fpga_counts_per_sample(fpga_counts_per_sample_),
    wt_cutoff(wt_cutoff_),
    target_gbps(target_gbps_),
    ibeam(_init_ivec(ibeam_,"beam")),
    ifreq_chunk(_init_ivec(ifreq_chunk_,"freq"))
{
    // Tons of argument checking.
    // The { nbeam, ibeam, nfreq_per_chunk, ifreq_chunk } args have already been checked in _init_ivec().

    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of nupfreq");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of fpga_counts_per_sample");

    if (nfreq_per_packet <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nfreq_per_packet > 0");
    if (nfreq_per_chunk % nfreq_per_packet)
	throw runtime_error("chime intensity_network_ostream constructor: expected nfreq_per_chunk to be a multiple of nfreq_per_packet");

    if (nt_per_chunk <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_chunk > 0");
    if (nt_per_packet <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_packet > 0");
    if (nt_per_chunk % nt_per_packet)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_chunk to be a multiple of nt_per_packet");

    if (wt_cutoff < 0.0)
	throw runtime_error("chime intensity_network_ostream constructor: expected wt_cutoff to be >= 0.0");
    if (target_gbps < 0.0)
	throw runtime_error("chime intensity_network_ostream constructor: expected target_gbps to be >= 0.0");

    this->npackets_per_chunk = (nfreq_per_chunk / nfreq_per_packet) * (nt_per_chunk / nt_per_packet);
    this->nbytes_per_packet = packet_size(nbeam, nfreq_per_packet, nupfreq, nt_per_packet);

    if (nbytes_per_packet > constants::max_output_udp_packet_size)
	throw runtime_error("chime intensity_network_ostream constructor: packet size is too large, you need to decrease nfreq_per_packet or nt_per_packet");

    xpthread_mutex_init(&this->state_lock);
    xpthread_cond_init(&this->cond_state_changed);

    int capacity = constants::output_ringbuf_capacity;
    string dropmsg = "warning: network write thread couldn't keep up with data, dropping packets";
    this->ringbuf = make_unique<udp_packet_ringbuf> (capacity, npackets_per_chunk, npackets_per_chunk * nbytes_per_packet, dropmsg);
    
    // Allocate encoding buffers
    int ndata = nbeam * nfreq_per_packet * nupfreq * nt_per_packet;
    this->tmp_intensity_vec.resize(ndata, 0.0);
    this->tmp_mask_vec.resize(ndata, 0.0);
    this->tmp_packet_list = ringbuf->allocate_packet_list();
}


intensity_network_ostream::~intensity_network_ostream()
{
    close(sockfd);
    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&state_lock);
}


void intensity_network_ostream::network_thread_startup()
{
    pthread_mutex_lock(&state_lock);
    network_thread_started = true;
    pthread_cond_broadcast(&cond_state_changed);
    pthread_mutex_unlock(&state_lock);
}


void intensity_network_ostream::wait_for_network_thread_startup()
{
    pthread_mutex_lock(&state_lock);
    
    while (!network_thread_started)
	pthread_cond_wait(&cond_state_changed, &state_lock);

    pthread_mutex_unlock(&state_lock);
}


// The 'intensity' and 'weights' arrays have shapes (nbeam, nfreq_per_chunk, nupfreq, nt_per_chunk)
void intensity_network_ostream::send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count,  bool is_blocking)
{
    if (fpga_count % fpga_counts_per_sample)
	throw runtime_error("intensity_network_ostream::send_chunk(): fpga count must be divisible by fpga_counts_per_sample");
    if (tmp_packet_list.curr_npackets > 0)
	throw runtime_error("intensity_network_ostream::send_chunk(): internal error: tmp_packet_list nonempty?!");

    // Pointers to contiguous arrays of shape (nbeam, nfreq_per_packet, nupfreq, nt_per_packet)
    float *tmp_intensity = &tmp_intensity_vec[0];
    float *tmp_mask = &tmp_mask_vec[0];

    // The number of packets per chunk is (nf_outer * nt_outer)
    int nf_outer = nfreq_per_chunk / nfreq_per_packet;
    int nt_outer = nt_per_chunk / nt_per_packet;

    // Outer loop over packets
    for (int if_outer = 0; if_outer < nf_outer; if_outer++) {
	for (int it_outer = 0; it_outer < nt_outer; it_outer++) {

	    // Copy input data into the { tmp_intensity, tmp_mask } arrays,
	    // in order to "de-stride" and apply the wt_cutoff.
	    
	    for (int ibeam = 0; ibeam < nbeam; ibeam++) {
		for (int if_inner = 0; if_inner < nfreq_per_packet; if_inner++) {
		    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
			// Row offset in 'tmp_intensity' and 'tmp_mask' arrays
			int idst = ibeam * nfreq_per_packet * nupfreq * nt_per_packet;
			idst += if_inner * nupfreq * nt_per_packet;
			idst += iupfreq * nt_per_packet;

			// Row offset in 'intensity' and 'weights' input arrays
			int isrc = ibeam * nfreq_per_chunk * nupfreq * stride;
			isrc += (if_outer * nfreq_per_packet + if_inner) * nupfreq * stride;
			isrc += iupfreq * stride;
			isrc += it_outer * nt_per_packet;

			// Operate on contiguous "row" of length nt_per_packet
			for (int it_inner = 0; it_inner < nt_per_packet; it_inner++) {
			    tmp_intensity[idst + it_inner] = intensity[isrc + it_inner];
			    tmp_mask[idst + it_inner] = (weights[isrc + it_inner] >= wt_cutoff) ? 1.0 : 0.0;
			}
		    }
		}
	    }

	    encode_packet(nbeam, nfreq_per_packet, nupfreq, nt_per_packet,
			  fpga_counts_per_sample, 
			  fpga_count + it_outer * nt_per_packet * fpga_counts_per_sample,
			  tmp_packet_list.data_end,                              // output buffer for encoding
			  &ibeam[0], &ifreq_chunk[if_outer * nfreq_per_packet],  // ibeam, ifreq arrays
			  tmp_intensity, tmp_mask);                              // input buffers for encoding

	    tmp_packet_list.add_packet(nbytes_per_packet);
	}
    }

    if (!ringbuf->producer_put_packet_list(tmp_packet_list, is_blocking)) 
	throw runtime_error("intensity_network_ostream::send_chunk() called after end_stream()");
}


void intensity_network_ostream::end_stream(bool join_network_thread)
{
    bool join_thread_after_releasing_lock = false;

    pthread_mutex_lock(&this->state_lock);

    if (join_network_thread && !network_thread_joined) {
	network_thread_joined = true;
	join_thread_after_releasing_lock = true;
    }

    pthread_mutex_unlock(&this->state_lock);

    this->ringbuf->end_stream();
    
    if (join_thread_after_releasing_lock)
	pthread_join(network_thread, NULL);
}


// static member function (de facto constructor)
auto intensity_network_ostream::make(const std::string &dstname, const std::vector<int> &ibeam_, 
				     const std::vector<int> &ifreq_chunk_, int nupfreq_, int nt_per_chunk_,
				     int nfreq_per_packet_, int nt_per_packet_, int fpga_counts_per_sample_,
				     float wt_cutoff_, double gbps_) -> std::shared_ptr<intensity_network_ostream>
{
    auto p = new intensity_network_ostream(dstname, ibeam_, ifreq_chunk_, nupfreq_, 
					   nt_per_chunk_, nfreq_per_packet_, nt_per_packet_, 
					   fpga_counts_per_sample_, wt_cutoff_, gbps_);
    
    shared_ptr<intensity_network_ostream> ret(p);
    xpthread_create(&ret->network_thread, network_thread_main, ret, "network write thread");
    ret->wait_for_network_thread_startup();

    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Network write thread


static void *network_thread_main(void *opaque_arg)
{
    auto stream = xpthread_get_arg<intensity_network_ostream> (opaque_arg, "network write thread");
    stream->network_thread_startup();

    try {
	network_thread_main2(stream.get());
    } catch (...) {
	stream->end_stream(false);   // "false" means "don't join threads" (would deadlock otherwise!)
	throw;
    }

    stream->end_stream(false);   // "false" has same meaning as above

    return NULL;
}


static void network_thread_main2(intensity_network_ostream *stream)
{
    int sockfd = stream->get_sockfd();
    double target_gbps = stream->get_target_gbps();
    udp_packet_list packet_list = stream->allocate_packet_list();

    // Used for throttling and displaying summary info at the end.
    struct timeval initial_time;
    int last_send_bytecount = 0;
    int64_t last_send_timestamp = 0;   // in usec
    int64_t npackets_sent = 0;
    int64_t nbytes_sent = 0;
    
    for (;;) {
	if (!stream->get_packet_list(packet_list))
	    break;   // end of stream reached (probably normal termination)
	
	for (int ipacket = 0; ipacket < packet_list.curr_npackets; ipacket++) {
	    const uint8_t *packet = packet_list.data_start + packet_list.packet_offsets[ipacket];
	    const int packet_nbytes = packet_list.packet_offsets[ipacket+1] - packet_list.packet_offsets[ipacket];
	    struct timeval curr_time = xgettimeofday();

	    if (npackets_sent == 0)
		initial_time = curr_time;

	    int64_t curr_timestamp = usec_between(initial_time, curr_time);

	    if ((target_gbps > 0.0) && (npackets_sent > 0)) {
		int64_t target_timestamp = last_send_timestamp + int64_t(8.0e-3 * last_send_bytecount / target_gbps);
		
		if (curr_timestamp < target_timestamp) {
		    xusleep(target_timestamp - curr_timestamp);
		    curr_timestamp = target_timestamp;
		}
	    }

	    // FIXME: sendmmsg() may improve performance here
	    ssize_t n = send(sockfd, packet, packet_nbytes, 0);
	    
	    if (n < 0)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() failed: ") + strerror(errno));
	    if (n != packet_nbytes)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() sent ") + to_string(n) + "/" + to_string(packet_nbytes) + " bytes?!");

	    last_send_bytecount = packet_nbytes;
	    last_send_timestamp = curr_timestamp;
	    nbytes_sent += packet_nbytes;
	    npackets_sent++;
	}
    }

    // FIXME temporary hack.  For testing, it is convenient to have a way of ending the stream. 
    // We make the rule that a packet with nbeam = nfreq = nupfreq = ntsamp = 0 means "end of stream".
    // Later this will be replaced by something better!  
    //
    // Since UDP doesn't guarantee delivery, we have no way to ensure that the end-of-stream packet 
    // reaches the other side, but we'll make a best effort by sending 10 packets separated by 0.1 sec.

    cerr << "ch_frb_io: network output thread sending end-of-stream packets\n";

    for (int ipacket = 0; ipacket < 10; ipacket++) {
	vector<uint8_t> packet(24, uint8_t(0));
	*((uint32_t *) &packet[0]) = uint32_t(1);  // protocol number

	ssize_t n = send(sockfd, &packet[0], packet.size(), 0);

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
	ss << ", gbps=" << (8.0e-3 * nbytes_sent / double(last_send_timestamp));
    if (target_gbps > 0.0)
	ss << ", target_gbps=" << target_gbps;
    ss << "\n";

    cerr << ss.str();
}


}  // namespace ch_frb_io
