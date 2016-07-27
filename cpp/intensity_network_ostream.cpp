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


// This is a conservative estimate for max UDP packet size without fragmentation, 
// on a 1Gbps ethernet link with jumbo frames enabled.
//
// FIXME a loose end here: is there a way to get the size at runtime?  In addition 
// to being more reliable, this would check to make sure jumbo frames are enabled.
//
// Note: setting this to a value >= 2^15 will cause problems with the current protocol,
// since we use a signed 16-bit packet length field.

static constexpr int max_packet_size = 8910;


// -------------------------------------------------------------------------------------------------
//
// chunk_exchanger helper class


struct chunk_exchanger {
    // Compile-time constants
    static constexpr int capacity = 4;
    static constexpr int unused_pool_size = capacity + 2;
    
    // Constant after constructor is called
    const int npackets_per_chunk;
    const int nbytes_per_packet;
    const int nbytes_per_chunk;

    // This socket is connected to the destination IP address and UDP port
    int sockfd;

    pthread_mutex_t mutex;
    pthread_cond_t cond_chunk_produced;
    pthread_cond_t cond_chunk_consumed;

    int curr_ipos;      // number of chunks consumed so far
    int curr_size;      // number of pending chunks in buffer
    uint8_t *curr_chunks[capacity];  // ring buffer

    int nunused;
    uint8_t *unused_chunks[unused_pool_size];

    bool endflag;

    chunk_exchanger(const std::string &dstname, int npackets_per_chunk_, int nbytes_per_packet_);
    ~chunk_exchanger();

    // Called by producer thread
    uint8_t *producer_put_chunk(uint8_t *chunk);
    
    // Called by consumer thread
    // Returns nullptr if stream has ended
    const uint8_t *consumer_get_chunk(const uint8_t *prev_chunk);

    void producer_end_stream();

    // Helper function called by constructor.  Returns file descriptor.
    static int make_socket_from_dstname(const std::string &dstname);
};


chunk_exchanger::chunk_exchanger(const std::string &dstname, int npackets_per_chunk_, int nbytes_per_packet_)
    : npackets_per_chunk(npackets_per_chunk_),
      nbytes_per_packet(nbytes_per_packet_),
      nbytes_per_chunk(npackets_per_chunk_ * nbytes_per_packet_),
      curr_ipos(0),
      curr_size(0),
      nunused(0),
      endflag(false)
{
    if (npackets_per_chunk < 1)
	throw runtime_error("chime intensity_network_ostream constructor: expected npackets_per_chunk >= 1");
    if (nbytes_per_packet < 1)
	throw runtime_error("chime intensity_network_ostream constructor: expected nbytes_per_packet >= 1");
    if (nbytes_per_packet > max_packet_size)
	throw runtime_error("chime intensity_network_ostream constructor: expected nbytes_per_packet <= " + to_string(max_packet_size));

    // Parse dstname and do DNS lookup if necessary.
    this->sockfd = make_socket_from_dstname(dstname);

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond_chunk_produced, NULL);
    pthread_cond_init(&cond_chunk_consumed, NULL);

    memset(curr_chunks, 0, sizeof(curr_chunks));
    memset(unused_chunks, 0, sizeof(unused_chunks));
}


chunk_exchanger::~chunk_exchanger()
{
    pthread_cond_destroy(&cond_chunk_produced);
    pthread_cond_destroy(&cond_chunk_consumed);
    pthread_mutex_destroy(&mutex);
    
    for (int i = 0; i < curr_size; i++) {
	int j = (curr_ipos + i) % capacity;
	free(curr_chunks[j]);
	curr_chunks[j] = nullptr;
    }

    for (int i = 0; i < nunused; i++) {
	free(unused_chunks[i]);
	unused_chunks[i] = nullptr;
    }

    close(sockfd);
    endflag = true;
}


uint8_t *chunk_exchanger::producer_put_chunk(uint8_t *old_chunk)
{
    pthread_mutex_lock(&mutex);

    if (endflag) {
	// Currently treated as an error, since the endflag is set by the producer thread
	pthread_mutex_unlock(&mutex);
	throw runtime_error("chime intensity_network_ostream: internal error: endflag is set in producer_put_chunk()?!");
    }

    if (old_chunk != NULL) {
	// Add old_chunk to current chunk list, blocking if necessary
	while (curr_size >= capacity)
	    pthread_cond_wait(&cond_chunk_consumed, &mutex);

	int i = (curr_ipos + curr_size) % capacity;
	curr_chunks[i] = old_chunk;
	curr_size++;

	pthread_cond_broadcast(&cond_chunk_produced);
    }

    uint8_t *ret = NULL;
    
    // Get new chunk from unused pool if possible...
    if (nunused > 0) {
	ret = unused_chunks[nunused];
	unused_chunks[nunused] = nullptr;
	nunused--;
    }

    pthread_mutex_unlock(&mutex);

    // ... if pool is empty, we call malloc() after releasing the lock
    if (!ret)
	ret = aligned_alloc<uint8_t> (nbytes_per_chunk);  // never returns NULL

    return ret;
}


const uint8_t *chunk_exchanger::consumer_get_chunk(const uint8_t *prev_chunk)
{
    const uint8_t *ret = nullptr;

    pthread_mutex_lock(&mutex);

    for (;;) {
	// return previous chunk to unused pool if possible
	if (prev_chunk && (nunused < unused_pool_size)) {
	    unused_chunks[nunused++] = const_cast<uint8_t *> (prev_chunk);
	    prev_chunk = nullptr;
	}

	// return pending chunk if it exists
	if (curr_size > 0) {
	    int i = curr_ipos % capacity;
	    ret = curr_chunks[i];
	    curr_chunks[i] = nullptr;
	    curr_ipos++;
	    curr_size--;

	    pthread_cond_broadcast(&cond_chunk_consumed);
	    break;
	}

	// end of stream (return NULL)
	if (endflag)
	    break;

	pthread_cond_wait(&cond_chunk_produced, &mutex);
    }
    
    pthread_mutex_unlock(&mutex);

    if (prev_chunk)
	free((void *) prev_chunk);
    
    return ret;
}


// static member function
int chunk_exchanger::make_socket_from_dstname(const string &dstname)
{
    // Parse dstname: expect string of the form HOSTNAME:PORT

    size_t i = dstname.find(":");

    if (i == std::string::npos)
	throw runtime_error("ch_frb_io: expected dstname='" + dstname + "' to be a colon-separated HOSTNAME:PORT string");
	
    string hostname = dstname.substr(0,i);
    string portstr = dstname.substr(i+1);
    uint16_t port = 0;

    try {
	port = lexical_cast<uint16_t> (portstr);
    }
    catch (...) {
	throw runtime_error("ch_frb_io: couldn't convert string '" + portstr + "' to 16-bit udp port number");
    }

    struct sockaddr_in saddr;
    memset(&saddr, 0, sizeof(saddr));
    saddr.sin_family = AF_INET;
    saddr.sin_port = htons(port);
    
    // FIXME need getaddrinfo() here
    int err = inet_pton(AF_INET, hostname.c_str(), &saddr.sin_addr);
    if (err == 0)
	throw runtime_error("ch_frb_io: couldn't resolve hostname '" + hostname + "' to an IP address: general parse error");
    if (err < 0)
	throw runtime_error("ch_frb_io: couldn't resolve hostname '" + hostname + "' to an IP address: " + strerror(errno) + "general parse error");

    int sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0)
	throw runtime_error(string("ch_frb_io: couldn't create udp socket: ") + strerror(errno));
    
    // Note: bind() not called, so source port number of outgoing packets will be arbitrarily assigned

    if (connect(sockfd, reinterpret_cast<struct sockaddr *> (&saddr), sizeof(saddr)) < 0) {
	close(sockfd);
	throw runtime_error("ch_frb_io: couldn't connect udp socket to dstname '" + dstname + "': " + strerror(errno));
    }

    return sockfd;
}


void chunk_exchanger::producer_end_stream()
{
    pthread_mutex_lock(&mutex);
    endflag = true;
    pthread_cond_broadcast(&cond_chunk_produced);
    pthread_mutex_unlock(&mutex);    
}


// -------------------------------------------------------------------------------------------------
//
// Network write thread


static void *network_thread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_thread_main()");

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<chunk_exchanger> *arg = (shared_ptr<chunk_exchanger> *) opaque_arg;
    shared_ptr<chunk_exchanger> exchanger = *arg;   // 'exchanger' is safe to use below
    delete arg;

    if (!exchanger)
	throw runtime_error("ch_frb_io: internal error: no chunk_exchanger passed to network_thread_main()");

    cerr << "ch_frb_io: network output thread starting\n";
    
    const int sockfd = exchanger->sockfd;
    const int npackets_per_chunk = exchanger->npackets_per_chunk;
    const int nbytes_per_packet = exchanger->nbytes_per_packet;
    const uint8_t *chunk = nullptr;
    ssize_t npackets_sent = 0;
    
    for (;;) {
	chunk = exchanger->consumer_get_chunk(chunk);
	
	// consumer_get_chunk() returns NULL if end of stream is reached (normal termination)
	if (!chunk)
	    break;
	
	// FIXME: sendmmsg() may improve performance here
	for (int ipacket = 0; ipacket < npackets_per_chunk; ipacket++) {
	    const uint8_t *packet = chunk + ipacket * nbytes_per_packet;
	    ssize_t n = send(sockfd, packet, nbytes_per_packet, 0);
	    
	    if (n < 0)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() failed: ") + strerror(errno));
	    if (n != nbytes_per_packet)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() sent ") + to_string(n) + "/" + to_string(nbytes_per_packet) + " bytes?!");

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

    cerr << ("ch_frb_io: network output thread exiting (" + to_string(npackets_sent) + " packets sent)\n");
    return NULL;
}



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
    float acc0 = 0.0;
    float acc1 = 0.0;
    float acc2 = 0.0;

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

    float mean = acc1/acc0;
    float var = acc2/acc0 - mean*mean;
    
    var = max(var, float(1.0e-10*mean*mean));

    float scale = sqrt(var) / 25.;
    float offset = -128.*scale + mean;   // 0x80 -> mean

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

    int nrows = nbeam * nfreq;
    int rowlen = nupfreq * ntsamp;

#if 0  // debug: print outputs
    for (int i = 0; i < nrows; i++) {
	cout << "XXX2\n";
	for (int j = 0; j < rowlen; j++)
	    cout << " " << intensity[i*rowlen+j];
	cout << "\n";
    }
#endif

    float *scale0 = (float *) (out + 24 + 2*nbeam + 2*nfreq);
    float *offset0 = (float *) (out + 24 + 2*nbeam + 2*nfreq + 4*nrows);
    uint8_t *data0 = out + 24 + 2*nbeam + 2*nfreq + 8*nrows;

    for (int irow = 0; irow < nrows; irow++)
	encode_packet_row(rowlen, scale0+irow, offset0+irow, data0 + irow*rowlen, intensity + irow*rowlen, mask + irow*rowlen);

#if 0  // debug
    int mask_count = 0;
    for (int i = 0; i < nrows*rowlen; i++) {
	if ((data0[i] == (uint8_t)0) || (data0[i] == (uint8_t)255))
	    mask_count++;
    }

    double mask_frac = (double)mask_count / (double)(nrows*rowlen);
    cout << ("XXX mask_frac = " + to_string(mask_frac) + "\n");
#endif

#if 0  // debug: print outputs
    for (int i = 0; i < nrows; i++) {
	cout << "XXX3\n";
	for (int j = 0; j < rowlen; j++)
	    cout << " " << (int)data0[i*rowlen+j];
	cout << "\n";
    }
#endif
}


// FIXME also in intensity_network_ostream.cpp
inline int packet_size(int nbeam, int nfreq, int nupfreq, int ntsamp)
{
    return 24 + 2*nbeam + 2*nfreq + 8*nbeam*nfreq + (nbeam * nfreq * nupfreq * ntsamp);
}


// -------------------------------------------------------------------------------------------------
//
// intensity_network_ostream


// Helper function for intensity_network_ostream constructor
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

    
intensity_network_ostream::intensity_network_ostream(const std::string &dstname, const std::vector<int> &ibeam_, 
						     const std::vector<int> &ifreq_chunk_, int nupfreq_, int nt_per_chunk_,
						     int nfreq_per_packet_, int nt_per_packet_, int fpga_counts_per_sample_,
						     float wt_cutoff_) :
    nbeam(ibeam_.size()),
    nupfreq(nupfreq_),
    nfreq_per_chunk(ifreq_chunk_.size()),
    nfreq_per_packet(nfreq_per_packet_),
    nt_per_chunk(nt_per_chunk_),
    nt_per_packet(nt_per_packet_),
    fpga_counts_per_sample(fpga_counts_per_sample_),
    wt_cutoff(wt_cutoff_),
    ibeam(_init_ivec(ibeam_,"beam")),
    ifreq_chunk(_init_ivec(ifreq_chunk_,"freq"))
{
    // Tons of argument checking.
    // The { nbeam, ibeam, nfreq_per_chunk, ifreq_chunk } args have already been checked in _init_ivec().

    if (nupfreq <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nupfreq > 0");
    if (nupfreq > 16)
	throw runtime_error("chime intensity_network_ostream constructor: expected nupfreq <= 16");

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

    if (fpga_counts_per_sample <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected fpga_counts_per_sample to be > 0");
    if (fpga_counts_per_sample >= 65536)
	throw runtime_error("chime intensity_network_ostream constructor: expected fpga_counts_per_sample to be < 2^16");

    if (wt_cutoff < 0.0)
	throw runtime_error("chime intensity_network_ostream constructor: expected wt_cutoff to be >= 0.0");

    this->nbytes_per_packet = packet_size(nbeam, nfreq_per_packet, nupfreq, nt_per_packet);

    if (nbytes_per_packet >= max_packet_size)
	throw runtime_error("chime intensity_network_ostream constructor: packet size is too large, you need to decrease nfreq_per_packet or nt_per_packet");

    // Allocate chunk_exchanger: a buffer for exchanging packets between the stream object and the network thread.
    // Note: the chunk_exchanger constructor resolves the dstname and creates a socket.

    int npackets_per_chunk = (nfreq_per_chunk / nfreq_per_packet) * (nt_per_chunk / nt_per_packet);
    this->exchanger = make_shared<chunk_exchanger> (dstname, npackets_per_chunk, nbytes_per_packet);

    // FIXME when code is reorganized this can go closer to the top of this function
    stringstream ss;
    ss << "ch_frb_io: constructing network_ostream (nbeams=" << nbeam << ",nfreq_per_chunk=" << nfreq_per_chunk
       << ",nt_per_chunk=" << nt_per_chunk << ",nfreq_per_packet=" << nfreq_per_packet << ",nt_per_packet=" << nt_per_packet
       << ",npackets_per_chunk=" << npackets_per_chunk << ")\n";

    string s = ss.str();
    cerr << s.c_str();
    
    // Allocate encoding buffers

    int ndata = nbeam * nfreq_per_packet * nupfreq * nt_per_packet;
    this->tmp_intensity_vec.resize(ndata, 0.0);
    this->tmp_mask_vec.resize(ndata, 0.0);

    this->packet_buf = exchanger->producer_put_chunk(nullptr);

    // Spawn network thread.
    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    
    shared_ptr<chunk_exchanger> *p = new shared_ptr<chunk_exchanger> (exchanger);
    int err = pthread_create(&this->network_thread, NULL, network_thread_main, p);

    if (err) {
	delete p;   // note: if pthread_create() succeeds, then slave thread will delete this pointer
	throw runtime_error(string("chime intensity_network_ostream constructor: couldn't create network thread: ") + strerror(errno));
    }

    this->network_thread_valid = true;
}


// The 'intensity' and 'weights' arrays have shapes (nbeam, nfreq_per_chunk, nupfreq, nt_per_chunk)
void intensity_network_ostream::send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count)
{
#if 0  // debug
    for (int i = 0; i < 64; i++) {
	cout << "XXX0";
	for (int j = 0; j < 64; j++)
	    cout << " " << intensity[i*stride+j];
	cout << "\n";
    }
#endif

    if (fpga_count % fpga_counts_per_sample)
	throw runtime_error("intensity_network_ostream::send_chunk(): fpga count must be divisible by fpga_counts_per_sample");

    // Pointers to contiguous arrays of shape (nbeam, nfreq_per_packet, nupfreq, nt_per_packet)
    float *tmp_intensity = &tmp_intensity_vec[0];
    float *tmp_mask = &tmp_mask_vec[0];

    // The number of packets per chunk is (nf_outer * nt_outer)
    int nf_outer = nfreq_per_chunk / nfreq_per_packet;
    int nt_outer = nt_per_chunk / nt_per_packet;

    // Outer loop over packets
    for (int if_outer = 0; if_outer < nf_outer; if_outer++) {
	for (int it_outer = 0; it_outer < nt_outer; it_outer++) {
	    int ipacket = if_outer * nt_outer + it_outer;

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

#if 0  // debug
	    for (int i = 0; i < 64; i++) {
		cout << "XXX1";
		for (int j = 0; j < 64; j++)
		    cout << " " << tmp_intensity[i*64+j];
		cout << "\n";
	    }
#endif

	    encode_packet(nbeam, nfreq_per_packet, nupfreq, nt_per_packet,
			  fpga_counts_per_sample, fpga_count,
			  packet_buf + ipacket * nbytes_per_packet,              // output buffer for encoding
			  &ibeam[0], &ifreq_chunk[if_outer * nfreq_per_packet],  // ibeam, ifreq arrays
			  tmp_intensity, tmp_mask);                              // input buffers for encoding
	}
    }

    // All packets have been encoded.  
    // Hand off packet buffer to network thread and get a new buffer.
    this->packet_buf = exchanger->producer_put_chunk(packet_buf);
}


void intensity_network_ostream::end_stream()
{
    exchanger->producer_end_stream();

    if (network_thread_valid) {
	pthread_join(network_thread, NULL);
	network_thread_valid = false;
    }
}


intensity_network_ostream::~intensity_network_ostream()
{
    if (packet_buf) {
	free(packet_buf);
	packet_buf = nullptr;
    }

    if (network_thread_valid) {
	pthread_cancel(network_thread);
	network_thread_valid = false;
    }
}


}  // namespace ch_frb_io
