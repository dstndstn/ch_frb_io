#ifndef _CH_FRB_IO_INTERNALS_HPP
#define _CH_FRB_IO_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include "ch_frb_io.hpp"


// This will compile to a "hint" for CPU branch prediction.  We use it mainly for error detection 
// in critical paths (e.g. packet parsing, where there are many possible ways a packet can be bad,
// leading to many unlikely branches in the code).  I found that it gives a few-percent speedup
// if used consistently, so it's not important, but seems worth doing anyway since it's so easy.

#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif


namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// struct intensity_packet: a lightweight struct representing one UDP packet.
// See L0_L1_packet.hpp for a more verbose description of the packet format.


struct intensity_packet {
    // "Header fields".   These 24 bytes should have the same ordering and byte count as the 
    // "on-wire" packet, since we use memcpy(..., 24) to initialize them from the raw packet data.
    uint32_t  protocol_version;
    int16_t   data_nbytes;
    uint16_t  fpga_counts_per_sample;
    uint64_t  fpga_count;
    uint16_t  nbeams;
    uint16_t  nfreq_coarse;
    uint16_t  nupfreq;
    uint16_t  ntsamp;

    // "Pointer" fields
    uint16_t  *beam_ids;          // 1D array of length nbeams
    uint16_t  *coarse_freq_ids;   // 1D array of length nfreq_coarse
    float     *scales;            // 2D array of shape (nbeam, nfreq_coarse)
    float     *offsets;           // 2D array of shape (nbeam, nfreq_coarse)
    uint8_t   *data;              // array of shape (nbeam, nfreq_coarse, nupfreq, ntsamp)


    static inline int packet_size(int nbeams, int nfreq_coarse, int nupfreq, int nt_per_packet)
    {
	int header_size = 24 + 2*nbeams + 2*nfreq_coarse + 8*nbeams*nfreq_coarse;
	int data_size = nbeams * nfreq_coarse * nupfreq * nt_per_packet;
	return header_size + data_size;
    }


    // Initializes a 'struct intensity_packet' from raw packet data.  The "pointer" fields of the
    // struct intensity_packet are initialized to pointers into the 'src' buffer, so the caller is
    // responsible for ensuring that this buffer doesn't get freed while the struct intensity_packet 
    // is in scope.
    //
    // Does a bunch of sanity checks and returns 'true' if packet is good, 'false' if bad.
    // (See extended comment in intensity_packet.cpp for a complete list of checks performed.)

    bool decode(const uint8_t *src, int src_nbytes);

    
    // Encodes a floating-point array of intensities into raw packet data, before sending packet.
    // The semantics of encode() aren't very intuitive, so we document them carefully here!
    //
    //    - Caller should initialize the "header" fields of the struct intensity packet.
    //
    //    - Caller should initialize the pointer fields 'beam_ids' and 'coarse_freq_ids' to 
    //      point to arrays of appropriate size.
    //
    //    - Caller should initialize the 'intensity' and 'weights' arrays to point to logical arrays 
    //      of shape (nbeams, nfreq_coarse, nupfreq, ntsamp).  This is the data that will be encoded 
    //      into the packet.  The stride arguments are defined so that the intensity array element with 
    //      logical indices (b,f,u,t) has memory location
    //
    //          intensity + b*beam_stride + (f*nupfreq+u)*freq_stride + t
    //
    //      and likewise for the weights.
    //
    //    - Since the binary packet format doesn't support weights but does support a boolean mask,
    //      encode() simply masksdata whose weight is below the 'wt_cutoff' argument.
    //
    //    - Caller must ensure that 'dst' points to a large enough buffer to encode the packet.
    //      For example, in intensity_network_ostream.cpp, we compute the buffer size ahead of
    //      time using intensity_packet::packet_size().
    //
    //    - Caller doesn't need to initialize the pointer fields 'scales', 'offsets', 'data'.
    //      These pointers are initialized in encode(), to point into the appropriate locations
    //      in the 'dst' buffer.  The actual scales and offsets are computed in encode() based
    //      on the mean and variance of the data.  The scales are chosen so that the intensities
    //      are masked if they deviate from the mean by approx 5 sigma.
    //
    // Returns size of the encoded packet in bytes.
    //
    // Caveat emptor: encode() doesn't do any argument checking at all, e.g. it's easy to segfault
    // by calling it wrong!

    int encode(uint8_t *dst, const float *intensity, const float *weights, int beam_stride, int freq_stride, float wt_cutoff);


    // Currently used only for debugging
    int find_coarse_freq_id(int id) const;
    bool contains_coarse_freq_id(int id) const;
};


// -------------------------------------------------------------------------------------------------
//
// udp_packet_list: a buffer containing opaque UDP packets.
// udp_packet_ringbuf: a thread-safe ring buffer for exchanging udp_packet_lists between threads.


struct udp_packet_list {
    // Capacity of buffer
    const int max_npackets;
    const int max_nbytes;  // summed over all packets

    // Current buffer size
    int curr_npackets = 0;
    int curr_nbytes = 0;   // summed over all packets
    bool is_full = false;

    // Packets are concatenated into the 'buf' array, and off_buf[i] stores the offset
    // of the i-th packet relative to the start of the buffer.  It's convenient to set
    // the sentinel value
    //   off_buf[curr_npackets] = curr_nbytes
    // so that the i-th packet always has size (off_buf[i+1] - off_buf[i])

    std::unique_ptr<uint8_t[]> buf;   // points to an array of length (max_nbytes + max_packet_size).
    std::unique_ptr<int[]> off_buf;   // points to an array of length (max_npackets + 1).

    // Bare pointers.
    uint8_t *data_start = nullptr;    // points to &buf[0]
    uint8_t *data_end = nullptr;      // points to &buf[curr_nbytes]
    int *packet_offsets = nullptr;    // points to &off_buf[0].  Note that packet_offsets[npackets] is always equal to 'nbytes'.

    udp_packet_list(int max_npackets, int max_nbytes);

    // Accessors (not range-checked)
    inline uint8_t *get_packet_data(int i)  { return data_start + packet_offsets[i]; }
    inline int get_packet_nbytes(int i)     { return packet_offsets[i+1] - packet_offsets[i]; }

    // To add a packet, we copy its data to the udp_packet_list::data_end pointer, then call add_packet()
    // to update the rest of the udp_packet_list fields consistently.
    void add_packet(int packet_nbytes);

    // Doesn't deallocate buffers or change the max_* fields, but sets the current packet count to zero.
    void reset();
};


// High-level comment: the get/put methods of udp_packet_ringbuf have been designed so that a
// fixed pool of udp_packet_lists is recycled throughout the lifetime of the ringbuf, rather than
// having buffers which are continually freed and allocated.  This is to avoid the page-faulting
// cost of Linux malloc.
//
// In contrast, the current implementation of assembled_chunk_ringbuf does lead to continually
// freed/allocated buffers.  It may improve performance to switch to an implementation which 
// is more similar to udp_packet_ringbuf!  (This loose end is noted in the README.)

struct udp_packet_ringbuf : noncopyable {
    // Specified at construction, used when new udp_packet_list objects are allocated
    const int ringbuf_capacity;
    const int max_npackets_per_list;
    const int max_nbytes_per_list;

    pthread_mutex_t lock;
    pthread_cond_t cond_packets_added;
    pthread_cond_t cond_packets_removed;
    bool stream_ended = false;

    int ringbuf_size = 0;
    int ringbuf_pos = 0;
    std::vector<std::unique_ptr<udp_packet_list> > ringbuf;

    udp_packet_ringbuf(int ringbuf_capacity, int max_npackets_per_list, int max_nbytes_per_list);
    ~udp_packet_ringbuf();
    
    // Note!  The pointer 'p' is _swapped_ with an empty udp_packet_list from the ring buffer.
    // In other words, when put_packet_list() returns, the argument 'p' points to an empty udp_packet_list.
    // Returns true on success, returns false if packets were dropped due to full ring buffer.
    // Throws an exception if called after end-of-stream.
    bool put_packet_list(std::unique_ptr<udp_packet_list> &p, bool is_blocking);

    // Note!  The pointer 'p' is _swapped_ with the udp_packet_list which is extracted from the ring buffer.
    // In other words, when get_packet_list() returns, the original udp_packet_list will be "recycled" (rather than freed).
    // Returns true on success (possibly after blocking), returns false if ring buffer is empty and stream has ended.
    bool get_packet_list(std::unique_ptr<udp_packet_list> &p);

    // Called by producer thread, when stream has ended.
    void end_stream();

    // No longer used, but I left it in anyway.
    bool is_alive();
};


// -------------------------------------------------------------------------------------------------
//
// assembled_chunk_ringbuf
//
// A thread-safe data structure for exchanging assembled_chunks between the assembler and processing threads.  
// It also manages the "active" assembled_chunks, which are being filled with data as new packets arrive.
// There is one assembled_chunk_ringbuf for each beam.


class assembled_chunk_ringbuf : noncopyable {
public:
    assembled_chunk_ringbuf(const intensity_network_stream::initializer &ini_params, int beam_id, int nupfreq,
			    int nt_per_packet, uint64_t fpga_counts_per_sample, uint64_t fpga_count0);

    ~assembled_chunk_ringbuf();

    // Called by assembler thread, to "assemble" an intensity_packet into the appropriate assembled_chunk.
    // The length-(intensity_network_stream::event_type::num_types) event_counts array is incremented 
    // in the appropriate indices.
    //
    // The packet must have nbeams==1.  If the actual packet stream has multiple beams, then each packet
    // is split into multiple packets (in a zero-copy way, see intensity_network_stream.cpp) which are
    // passed to different assembled_chunk_ringbufs.
    //
    // This routine is nonblocking.  If it would need to block, then it drops data somewhere and
    // increments the appropriate event_count.  (Depending which flags are set in ini_params, it may
    // also print a warning or throw an exception.)

    void put_unassembled_packet(const intensity_packet &packet, int64_t *event_counts);
    
    // Called when the assembler thread exits.  
    // Moves any remaining active chunks into the ring buffer and sets 'doneflag'.
    void end_stream(int64_t *event_counts);

    // Called by processing threads, via intensity_network_stream::get_assembled_chunk().
    // Returns the next assembled_chunk from the ring buffer, blocking if necessary to wait for data.
    // If the ring buffer is empty and end_stream() has been called, it returns an empty pointer
    // to indicate end-of-stream.
    std::shared_ptr<assembled_chunk> get_assembled_chunk();


protected:
    const intensity_network_stream::initializer ini_params;

    const int beam_id;
    const int nupfreq;
    const int nt_per_packet;
    const uint64_t fpga_counts_per_sample;

    // Helper function: adds assembled chunk to the ring buffer
    void _put_assembled_chunk(const std::shared_ptr<assembled_chunk> &chunk, int64_t *event_counts);

    // Helper function: allocates new assembled chunk
    std::shared_ptr<assembled_chunk> _make_assembled_chunk(uint64_t ichunk);

    // The "active" chunks are in the process of being filled with data as packets arrive.
    // Currently we take the active window to be two assembled_chunks long, but this could be generalized.
    // When an active chunk is finished, it is added to the ring buffer.
    // Note: the active_chunk pointers are not protected by a lock, but are only accessed by the assembler thread.
    // Note: active_chunk0->ichunk is always equal to (assembled_ringbuf_pos + assembled_ringbuf_size).
    std::shared_ptr<assembled_chunk> active_chunk0;
    std::shared_ptr<assembled_chunk> active_chunk1;

    // Not sure if this really affects bottom-line performance, but thought it would be a good idea
    // to ensure that the "assembler-only" and "shared" fields were on different cache lines.
    char pad[constants::cache_line_size];

    // All fields below are protected by the lock
    pthread_mutex_t lock;

    // Processing thread waits here if the ring buffer is empty.
    pthread_cond_t cond_assembled_chunks_added;

    std::shared_ptr<assembled_chunk> assembled_ringbuf[constants::assembled_ringbuf_capacity];
    int assembled_ringbuf_pos = 0;
    int assembled_ringbuf_size = 0;
    bool doneflag = false;
};


// -------------------------------------------------------------------------------------------------
//
// Miscelleanous inlines


inline bool is_power_of_two(int n)
{
    if (n <= 0)
	throw std::runtime_error("ch_frb_io: internal error: is_power_of_two() received argument <= 0");
    return (n & (n-1)) == 0;
}

inline int round_down_to_power_of_two(int n)
{
    if (n <= 0)
	throw std::runtime_error("ch_frb_io: internal error: is_power_of_two() received argument <= 0");
    return 1 << (int)log2(n+0.5);
}

inline int randint(std::mt19937 &rng, int lo, int hi)
{
    return std::uniform_int_distribution<>(lo,hi-1)(rng);   // note hi-1 here!
}

inline double uniform_rand(std::mt19937 &rng)
{
    return std::uniform_real_distribution<>()(rng);
}

inline double uniform_rand(std::mt19937 &rng, double lo, double hi)
{
    return lo + (hi-lo) * uniform_rand(rng);
}

template<typename T> inline void uniform_rand(std::mt19937 &rng, T *p, int n)
{
    for (int i = 0; i < n; i++)
	p[i] = uniform_rand(rng);
}

inline bool file_exists(const std::string &filename)
{
    struct stat s;

    int err = stat(filename.c_str(), &s);
    if (err >= 0)
        return true;
    if (errno == ENOENT)
        return false;

    throw std::runtime_error(filename + ": " + strerror(errno));
}

template<typename T> inline T prod(const std::vector<T> &v)
{
    T ret = (T)1;
    for (unsigned int i = 0; i < v.size(); i++)
	ret *= v[i];
    return ret;
}

// returns string representation of a vector
template<typename T> inline std::string vstr(const T *buf, int n, int stride=1)
{
    std::stringstream ss;
    ss << "[";
    for (int i = 0; i < n; i++)
	ss << " " << buf[i*stride];
    ss << " ]";
    return ss.str();
}

template<typename T> inline std::string vstr(const std::vector<T> &buf)
{
    return vstr(&buf[0], buf.size());
}


template<typename T>
inline T *aligned_alloc(size_t nelts)
{
    if (nelts == 0)
	return NULL;

    // align to 64-byte cache lines
    void *p = NULL;
    if (posix_memalign(&p, 64, nelts * sizeof(T)) != 0)
	throw std::runtime_error("couldn't allocate memory");

    memset(p, 0, nelts * sizeof(T));
    return reinterpret_cast<T *> (p);
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&& ...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


inline struct timeval xgettimeofday()
{
    struct timeval tv;

    int err = gettimeofday(&tv, NULL);
    if (err)
	throw std::runtime_error("gettimeofday failed");

    return tv;
}

inline int64_t usec_between(const struct timeval &tv1, const struct timeval &tv2)
{
    return 1000000 * int64_t(tv2.tv_sec - tv1.tv_sec) + int64_t(tv2.tv_usec - tv1.tv_usec);
}

inline void xusleep(useconds_t usec)
{
    int err = usleep(usec);
    if (err)
	throw std::runtime_error("usleep failed");
}

}  // namespace ch_frb_io

#endif // _CH_FRB_IO_INTERNALS_HPP
