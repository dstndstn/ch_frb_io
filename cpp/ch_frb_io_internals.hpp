#ifndef _CH_FRB_IO_INTERNALS_HPP
#define _CH_FRB_IO_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <cmath>
#include <random>
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


inline int randint(std::mt19937 &rng, int lo, int hi)
{
    return std::uniform_int_distribution<>(lo,hi+1)(rng);   // note hi+1 here!
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
template<typename T> inline std::string vstr(const T *buf, int n)
{
    std::stringstream ss;
    ss << "[";
    for (int i = 0; i < n; i++)
	ss << " " << buf[i];
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

inline int64_t usec_between(struct timeval &tv1, struct timeval &tv2)
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

template<typename T> inline void xpthread_create(pthread_t *thread, void *(*thread_main)(void *), const std::shared_ptr<T> &arg, const std::string &thread_name)
{
    // To pass a shared_ptr to a new pthread, we use a bare pointer to a shared_ptr.
    std::shared_ptr<T> *p = new std::shared_ptr<T> (arg);

    // If pthread_create() succeeds, then spawned thread is responsible for 'delete p'
    int err = pthread_create(thread, NULL, thread_main, p);

    if (err) {
	delete p;
	throw std::runtime_error("couldn't create " + thread_name + ": " + strerror(errno));
    }
}

template<typename T> inline std::shared_ptr<T> xpthread_get_arg(void *opaque_arg, const std::string &thread_name)
{
    if (!opaque_arg)
	throw std::runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to " + thread_name);

    std::shared_ptr<T> *arg = (std::shared_ptr<T> *) opaque_arg;
    std::shared_ptr<T> ret = *arg;
    delete arg;

    if (!ret)
	throw std::runtime_error("ch_frb_io: internal error: empty pointer passed to " + thread_name);

    return ret;
}


// Utility routine: converts a string to type T (only a few T's are defined; see lexical_cast.cpp)
template<typename T> extern T lexical_cast(const std::string &x);

// Unit test
extern void test_lexical_cast();


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_INTERNALS_HPP
