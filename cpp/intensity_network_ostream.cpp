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
    const double throughput_target;

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

    chunk_exchanger(const std::string &dstname, int npackets_per_chunk_, int nbytes_per_packet_, double throughput_target_);
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


chunk_exchanger::chunk_exchanger(const std::string &dstname, int npackets_per_chunk_, int nbytes_per_packet_, double throughput_target_)
    : npackets_per_chunk(npackets_per_chunk_),
      nbytes_per_packet(nbytes_per_packet_),
      nbytes_per_chunk(npackets_per_chunk_ * nbytes_per_packet_),
      throughput_target(throughput_target_),
      curr_ipos(0),
      curr_size(0),
      nunused(0),
      endflag(false)
{
    if (npackets_per_chunk < 1)
	throw runtime_error("chime intensity_network_ostream constructor: expected npackets_per_chunk >= 1");
    if (nbytes_per_packet < 1)
	throw runtime_error("chime intensity_network_ostream constructor: expected nbytes_per_packet >= 1");
    if (throughput_target < 0.01)
	throw runtime_error("chime intensity_network_ostream constructor: expected throughput_target >= 0.01");

    //
    // 8910 bytes is a conservative estimate for max UDP payload on a 1 Gbps ethernet link with jumbo frames enabled.
    //
    // FIXME a loose end here: is there a way to get the max payload size at runtime?  In addition to being more
    // reliable, this would check to make sure jumbo frames are enabled.
    //
    if (nbytes_per_packet > 8910)
	throw runtime_error("chime intensity_network_ostream constructor: expected nbytes_per_packet < 8900");

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
	ret = aligned_alloc<uint8_t> (nbytes_per_chunk);

    return ret;
}


const uint8_t *chunk_exchanger::consumer_get_chunk(const uint8_t *prev_chunk)
{
    const uint8_t *ret = nullptr;

    pthread_mutex_lock(&mutex);

    for (;;) {
	if (prev_chunk && (nunused < unused_pool_size)) {
	    unused_chunks[nunused++] = const_cast<uint8_t *> (prev_chunk);
	    prev_chunk = nullptr;
	}

	if (endflag)
	    break;

	if (curr_size > 0) {
	    int i = curr_ipos % capacity;
	    ret = curr_chunks[i];
	    curr_chunks[i] = nullptr;
	    curr_ipos++;
	    curr_size--;

	    pthread_cond_broadcast(&cond_chunk_consumed);
	    break;
	}

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
    string hostname;
    uint16_t port = 0;

    // Parse dstname: expect string of the form HOSTNAME:PORT

    try {
	size_t i = dstname.find(":");

	if (i == std::string::npos)
	    throw runtime_error("ch_frb_io: couldn't parse dstname='" + dstname + "'");
	
	hostname = dstname.substr(0,i);
	port = lexical_cast<uint16_t> (dstname.substr(i+1));
    }
    catch (...) {
	throw runtime_error("ch_frb_io: couldn't parse dstname='" + dstname + "'");
    }

    struct sockaddr_in saddr;
    memset(&saddr, 0, sizeof(saddr));
    saddr.sin_family = AF_INET;
    saddr.sin_port = htons(port);
    
    int err = inet_pton(AF_INET, hostname.c_str(), &saddr.sin_addr);
    if (err == 0)
	throw runtime_error("ch_frb_io: couldn't parse dstname='" + dstname + "'");
    if (err < 0)
	throw runtime_error("ch_frb_io: couldn't parse dstname='" + dstname + "': " + strerror(errno));

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
    pthread_mutex_unlock(&mutex);    
}


// -------------------------------------------------------------------------------------------------
//
// Network write thread
//
// FIXME implement throughput_target


static void *network_thread_main(void *opaque_arg)
{
    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr
    shared_ptr<chunk_exchanger> *arg = (shared_ptr<chunk_exchanger> *) opaque_arg;
    shared_ptr<chunk_exchanger> exchanger = *arg;   // 'exchanger' is safe to use below
    delete arg;
    
    const int sockfd = exchanger->sockfd;
    const int npackets_per_chunk = exchanger->npackets_per_chunk;
    const int nbytes_per_packet = exchanger->nbytes_per_packet;
    const double throughput_target = exchanger->throughput_target;

    const uint8_t *chunk = nullptr;
    
    for (;;) {
	chunk = exchanger->consumer_get_chunk(chunk);
	if (!chunk)
	    return NULL;
       
	// FIXME: sendmmsg() may improve performance here
	for (int ipacket = 0; ipacket < npackets_per_chunk; ipacket++) {
	    const uint8_t *packet = chunk + ipacket * nbytes_per_packet;
	    ssize_t n = send(sockfd, packet, nbytes_per_packet, 0);

	    if (n < 0)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() failed:") + strerror(errno));
	    if (n != nbytes_per_packet)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() sent ") + to_string(n) + "/" + to_string(nbytes_per_packet) + " bytes?!");
	}
    }
}


}  // namespace ch_frb_io
