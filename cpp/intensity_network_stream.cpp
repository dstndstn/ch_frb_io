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
// intensity_packet_stream


bool intensity_packet_stream::process_packet(int packet_nbytes, const uint8_t *in)
{
    if (packet_nbytes < 32)
	return false;
    
    uint32_t protocol_version = *((uint32_t *) in);
    if (protocol_version != 1)
	return false;

    int data_nbytes = *((int16_t *) (in+4));
    int nbeam = *((uint16_t *) (in+16));
    int nfreq = *((uint16_t *) (in+18));
    int nupfreq = *((uint16_t *) (in+20));
    int ntsamp = *((uint16_t *) (in+22));

    // The following way of writing the comparisons (using uint64_t) guarantees no overflow
    uint64_t n4 = uint64_t(nbeam) * uint64_t(nfreq) * uint64_t(nupfreq) * uint64_t(ntsamp);
    if (!n4 || (n4 > 10000))
	return false;
    if (data_nbytes != int(n4))
	return false;
    if (packet_nbytes != packet_size(nbeam,nfreq,nupfreq,ntsamp))
	return false;

    cerr << "!!! good packet: nbeam=" << nbeam << ", nfreq=" << nfreq << ", nupfreq=" << nupfreq << ", ntsamp=" << ntsamp << " !!!\n";
    return true;
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

    cerr << "waiting for packets...\n";
    
    for (;;) {
	ssize_t bytes_read = read(sock_fd, (char *) &packet_buf[0], max_packet_size);

	if (bytes_read >= max_packet_size) {
	    cerr << "!!! bad packet !!!\n";
	    continue;
	}

	if (!stream->process_packet(bytes_read, &packet_buf[0]))
	    cerr << "!!! bad packet !!!\n";
    }
}


void spawn_network_thread(int udp_port, const shared_ptr<intensity_packet_stream> &stream)
{
    if ((udp_port <= 0) || (udp_port >= 65536))
	throw runtime_error("ch_frb_io: spawn_network_thread(): bad udp port " + to_string(udp_port));

    shared_ptr<network_thread_context> context = make_shared<network_thread_context> (udp_port, stream);

    // To hand off the context to the network thread, we use a bare pointer to a shared pointer
    shared_ptr<network_thread_context> *p = new shared_ptr<network_thread_context> (context);
    pthread_t network_thread;

    int err = pthread_create(&network_thread, NULL, network_thread_main, p);
    
    if (err) {
	delete p;   // thread didn't start, caller is responsible for deleting p
	throw runtime_error("ch_frb_io: spawn_network_thread(): pthread_create failed");
    }
}


}  // namespace ch_frb_io
