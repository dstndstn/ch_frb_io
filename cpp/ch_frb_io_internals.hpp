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


// Branch predictor hint
#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif


namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------


struct intensity_packet {
    // "Header" fields
    uint32_t  protocol_version;
    int16_t   data_nbytes;
    uint16_t  fpga_counts_per_sample;
    uint64_t  fpga_count;
    uint16_t  nbeams;
    uint16_t  nfreq_coarse;
    uint16_t  nupfreq;
    uint16_t  ntsamp;

    // "Data" fields
    uint16_t  *beam_ids;   // 1D array of length nbeams
    uint16_t  *freq_ids;   // 1D array of length nfreq_coarse
    float     *scales;     // 2D array of shape (nbeam, nfreq_coarse)
    float     *offsets;    // 2D array of shape (nbeam, nfreq_coarse)
    uint8_t   *data;       // array of shape (nbeam, nfreq_coarse, nupfreq, ntsamp)

    // FIXME rethink the member function names below, since they're not very intuitive

    // Returns true if packet is good, false if bad
    bool read(const uint8_t *src, int src_nbytes);

    // Returns packet_nbytes (not data_nbytes)
    int write(uint8_t *dst) const;

    // The semantics of encode() aren't particularly intuitive, so we document it carefully here!
    //
    // It is assumed that the caller has initialized the "header" fields and the data pointers 'beam_ids', 'freq_ids'.
    // The other three data pointers (scales, offsets, data) are not assumed initialized, and encode() will initialize them
    // to point to subregions of the 'dst' array.
    //
    // The 'dst' array should point to an allocated but uninitialized memory region which will be filled by encode_packet().
    // Note that encode_packet() doesn't check for overflows (or do any arugment checking at all!)
    //
    // The 'intensity' and 'weights' pointers should point to arrays of logical shape (nbeams, nfreq_coarse, nupfreq, ntsamp).
    //
    // The stride arguments are defined so that the intensity array element with logical indices (b,f,u,t) is stored at
    // memory location
    //
    //    intensity + b*beam_stride + (f*nupfreq+u)*freq_stride + t
    //
    // and likewise for the weights.

    void encode(uint8_t *dst, const float *intensity, const float *weights, int beam_stride, int freq_stride, float wt_cutoff);

    // Currently used only for debugging
    int find_freq_id(int freq_id) const;
    bool contains_freq_id(int freq_id) const;
};


// FIXME make these static member functions of intensity_packet?
inline int header_size(int nbeams, int nfreq_coarse)
{
    return 24 + 2*nbeams + 2*nfreq_coarse + 8*nbeams*nfreq_coarse;
}

inline int packet_size(int nbeams, int nfreq_coarse, int nupfreq, int nt_per_packet)
{
    return header_size(nbeams, nfreq_coarse) + (nbeams * nfreq_coarse * nupfreq * nt_per_packet);
}


// -------------------------------------------------------------------------------------------------


struct udp_packet_ringbuf : noncopyable {
    pthread_mutex_t lock;
    pthread_cond_t cond_packets_added;
    pthread_cond_t cond_packets_removed;
    bool stream_ended = false;
    bool drops_allowed = true;

    const int ringbuf_capacity;
    int ringbuf_size = 0;
    int ringbuf_pos = 0;
    std::vector<udp_packet_list> ringbuf;

    // Specified at construction, used when new udp_packet_list objects are allocated
    const int max_npackets_per_list = 0;
    const int max_nbytes_per_list = 0;
    
    // This message is printed to stderr whenever packets are dropped (if empty string, no message will be printed)
    std::string dropmsg;

    // If 'drops_allowed' is false, the process will crash if the ring buffer overfills.
    // This sometimes makes sense during testing but probably not otherwise.
    udp_packet_ringbuf(int ringbuf_capacity, int max_npackets_per_list, int max_nbytes_per_list, 
		       const std::string &dropmsg = std::string(), bool drops_allowed = true);

    ~udp_packet_ringbuf();
    
    //
    // Important note!  These routines _swap_ their udp_packet_list argument with a packet_list in the ringbuf.
    //
    // I.e., producer_put_packet_list() is called with a full packet list, and swaps it for an empty packet list, which the
    // producer thread can fill with packets.  Similary, consumer_get_packet_list() is called with a "junk" packet list (which
    // may or may not be empty), and swaps it for a full packet list.
    //
    // The producer_put_packet_list() routine is nonblocking; if the ring buffer is full then it prints an error message 
    // ("dropmsg") and discards the packets.  Both routines return false if stream has ended, true otherwise.
    // 
    bool producer_put_packet_list(udp_packet_list &l, bool is_blocking);
    bool consumer_get_packet_list(udp_packet_list &l);
    
    // Can be called by either producer or consumer thread
    void end_stream();
    bool is_alive();
};


// -------------------------------------------------------------------------------------------------


class assembled_chunk_ringbuf : noncopyable {
public:
    // When the assembler is constructed, the fp_ fields must be initialized.
    assembled_chunk_ringbuf(const intensity_network_stream &s, int assembler_ix);
    ~assembled_chunk_ringbuf();

    // Called by assembler thread
    bool put_unassembled_packet(const intensity_packet &packet);
    void end_stream();   // called when assembler thread exits

    // Called by "processing" threads, via intensity_network_stream::get_assembled_chunk().
    std::shared_ptr<assembled_chunk> get_assembled_chunk();


protected:
    // Initialized at construction
    const intensity_network_stream::initializer _initializer;

    int beam_id = 0;
    int nupfreq = 0;
    int nt_per_packet = 0;
    int fpga_counts_per_sample = 0;

    std::shared_ptr<assembled_chunk> _make_assembled_chunk(uint64_t chunk_t0);
    void _put_assembled_chunk(const std::shared_ptr<assembled_chunk> &chunk);

    // This data is not protected by the lock, but is only accessed by the assembler thread.
    std::shared_ptr<assembled_chunk> active_chunk0;
    std::shared_ptr<assembled_chunk> active_chunk1;

    // All state below is protected by the lock
    pthread_mutex_t lock;
    pthread_cond_t cond_assembled_chunks_added;

    std::shared_ptr<assembled_chunk> assembled_ringbuf[constants::assembled_ringbuf_capacity];
    int assembled_ringbuf_pos = 0;
    int assembled_ringbuf_size = 0;
    bool doneflag = false;
};


// -------------------------------------------------------------------------------------------------


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

inline void xpthread_mutex_init(pthread_mutex_t *lock)
{
    if (pthread_mutex_init(lock, NULL) != 0)
        throw std::runtime_error("pthread_mutex_init() failed?!");    	    
}

inline void xpthread_cond_init(pthread_cond_t *cond)
{
    if (pthread_cond_init(cond, NULL) != 0)
        throw std::runtime_error("pthread_cond_init() failed?!");    	    
}


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_INTERNALS_HPP
