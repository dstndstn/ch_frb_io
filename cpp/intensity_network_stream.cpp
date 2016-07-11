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


struct network_thread_context {
    uint16_t udp_port;

    network_thread_context(uint16_t udp_port_) :
	udp_port(udp_port_)
    { }
};


static void *network_thread_main(void *opaque_arg)
{
    static const int socket_bufsize = 2 << 26;   // 64 MB

    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_thread_main()");

    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    shared_ptr<network_thread_context> *arg = (shared_ptr<network_thread_context> *) opaque_arg;
    shared_ptr<network_thread_context> context = *arg;   // 'context' is safe to use below
    delete arg;

    if (!context)
	throw runtime_error("ch_frb_io: internal error: no network_thread_context passed to network_thread_main()");
    
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

    int max_packet_size = 16384;
    vector<char> packet_buf(max_packet_size);

    cerr << "waiting for packets...\n";
    
    for (;;) {
	ssize_t bytes_read = read(sock_fd, &packet_buf[0], max_packet_size);
	cerr << "!!! received packet: " << bytes_read << " bytes !!!\n";
    }
}


void spawn_network_istream(int udp_port)
{
    if ((udp_port <= 0) || (udp_port >= 65536))
	throw runtime_error("spawn_network_istream(): bad udp port " + to_string(udp_port));

    shared_ptr<network_thread_context> context = make_shared<network_thread_context> (udp_port);

    // To hand off the context to the network thread, we use a bare pointer to a shared pointer
    shared_ptr<network_thread_context> *p = new shared_ptr<network_thread_context> (context);
    pthread_t network_thread;

    int err = pthread_create(&network_thread, NULL, network_thread_main, p);
    
    if (err) {
	delete p;   // thread didn't start, caller is responsible for deleting p
	throw runtime_error("spawn_network_istream: pthread_create() failed");
    }
}


}  // namespace ch_frb_io
