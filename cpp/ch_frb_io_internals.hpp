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

    //
    // Returns true if packet is good, false if bad
    //
    // Includes the following checks:
    //   - protocol version == 1
    //   - dimensions (nbeams, nfreq_coarse, nupfreq, ntsamp) are not large enough to lead to integer overflows
    //   - packet and data byte counts are correct
    //
    // Does not check the following:
    //   - any checks on beam ids or coarse_freq_ids
    //   - ntsamp is a power of two
    //   - nbeams, nfreq_coarse, nupfreq, ntsamp, fpga_counts_per_sample are all > 0
    //   - fpga_count is a multiple of (fpga_counts_per_sample * ntsamp)
    //
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


inline int header_size(int nbeams, int nfreq_coarse)
{
    return 24 + 2*nbeams + 2*nfreq_coarse + 8*nbeams*nfreq_coarse;
}

inline int packet_size(int nbeams, int nfreq_coarse, int nupfreq, int nt_per_packet)
{
    return header_size(nbeams, nfreq_coarse) + (nbeams * nfreq_coarse * nupfreq * nt_per_packet);
}


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
