#ifndef _CH_FRB_IO_INTERNALS_HPP
#define _CH_FRB_IO_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <hdf5.h>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include "ch_frb_io.hpp"

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif


inline double uniform_rand()
{
    return (rand() + 0.5) / (RAND_MAX + 1.0);
}

inline int randint(int lo, int hi)
{
    int ret = lo + (int)((hi-lo)*uniform_rand());
    ret = std::max(ret, lo);    // should be redundant
    ret = std::min(ret, hi-1);  // should be redundant
    return ret;
}

template<typename T> inline void uniform_rand(T *p, int n)
{
    for (int i = 0; i < n; i++)
	p[i] = uniform_rand();
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


// -------------------------------------------------------------------------------------------------
//
// HDF5 wrappers (these are generally useful since the libhdf5 C/C++ api is so clunky)


template<typename T> inline hid_t hdf5_type();

// Reference: https://www.hdfgroup.org/HDF5/doc/H5.user/Datatypes.html
template<> inline hid_t hdf5_type<int>()            { return H5T_NATIVE_INT; }
template<> inline hid_t hdf5_type<float>()          { return H5T_NATIVE_FLOAT; }
template<> inline hid_t hdf5_type<double>()         { return H5T_NATIVE_DOUBLE; }
template<> inline hid_t hdf5_type<unsigned char>()  { return H5T_NATIVE_UCHAR; }


struct hdf5_file : noncopyable {
    std::string filename;
    hid_t file_id;

    // If write=false, the file is opened read-only, and an exception is thrown if it doesn't exist.
    // If write=true, the file is opened for writing.  If the file already exists, it will either be clobbered
    // or an exception will be thrown, depending on the value of 'clobber'.
    hdf5_file(const std::string &filename, bool write=false, bool clobber=true);
    ~hdf5_file();
};


struct hdf5_group : noncopyable {
    std::string filename;
    std::string group_name;
    hid_t group_id;

    // If create=true, the group will be created if it doesn't exist.
    hdf5_group(const hdf5_file &f, const std::string &group_name, bool create=false);
    ~hdf5_group();

    bool has_attribute(const std::string &attr_name) const;
    bool has_dataset(const std::string &dataset_name) const;

    void get_attribute_shape(const std::string &attr_name, std::vector<hsize_t> &shape) const;
    void get_dataset_shape(const std::string &attr_name, std::vector<hsize_t> &shape) const;
    
    // For a string-valued dataset, 'slen' is the string length, plus one byte for the null terminator
    ssize_t get_dataset_slen(const std::string &dataset_name) const;

    // Read scalar attribute
    template<typename T> T read_attribute(const std::string &attr_name) const
    {
	T ret;
	this->_read_attribute(attr_name, hdf5_type<T>(), reinterpret_cast<void *> (&ret), std::vector<hsize_t>());
	return ret;
    }

    // Write scalar attribute
    template<typename T> void write_attribute(const std::string &attr_name, const T &x)
    {
	this->_write_attribute(attr_name, hdf5_type<T>(), reinterpret_cast<const void *> (&x), std::vector<hsize_t>());
    }

    // Write 1D vector attribute
    template<typename T> void write_attribute(const std::string &attr_name, const std::vector<T> &x)
    {
	std::vector<hsize_t> shape(1, x.size());
	this->_write_attribute(attr_name, hdf5_type<T>(), reinterpret_cast<const void *> (&x[0]), shape);
    }

    // Read multidimensional dataset
    template<typename T> void read_dataset(const std::string &dataset_name, T *out, const std::vector<hsize_t> &expected_shape) const
    {
	this->_read_dataset(dataset_name, hdf5_type<T>(), reinterpret_cast<void *> (out), expected_shape);
    }

    // Read multidimensional string-valued dataset.
    // The 'slen' arg should be the string length, plus one for the null terminator (as returned by hdf5_group::get_dataset_slen()).
    // The 'out' arg should point to an array of length prod(expected_shape) * slen.
    void read_string_dataset(const std::string &dataset_name, char *out, const std::vector<hsize_t> &expected_shape, ssize_t slen) const;
    
    // Write multidimensional dataset
    template<typename T> void write_dataset(const std::string &dataset_name, const T *data, const std::vector<hsize_t> &shape)
    {
	this->_write_dataset(dataset_name, hdf5_type<T>(), reinterpret_cast<const void *> (data), shape);
    }

    // Helpers
    void _get_attribute_shape(const std::string &attr_name, hid_t attr_id, std::vector<hsize_t> &shape) const;
    void _read_attribute(const std::string &attr_name, hid_t hdf5_type, void *out, const std::vector<hsize_t> &expected_shape) const;
    void _write_attribute(const std::string &attr_name, hid_t hdf5_type, const void *data, const std::vector<hsize_t> &shape);
    void _get_dataset_shape(const std::string &dataset_name, hid_t dataset_id, std::vector<hsize_t> &shape) const;
    void _read_dataset(const std::string &dataset_name, hid_t hdf5_type, void *out, const std::vector<hsize_t> &expected_shape) const;
    void _write_dataset(const std::string &dataset_name, hid_t hdf5_type, const void *data, const std::vector<hsize_t> &shape);
};


// This class isn't intended to be used directly; use the wrapper hdf5_extendable_dataset<T> below
struct _hdf5_extendable_dataset : noncopyable {
    std::string filename;
    std::string group_name;
    std::string dataset_name;
    std::string full_name;
    std::vector<hsize_t> curr_shape;
    int axis;

    hid_t type;
    hid_t dataset_id;

    _hdf5_extendable_dataset(const hdf5_group &g, const std::string &dataset_name, 
			     const std::vector<hsize_t> &chunk_shape, int axis, hid_t type, int bitshuffle);

    ~_hdf5_extendable_dataset();

    void write(const void *data, const std::vector<hsize_t> &shape);
};


template<typename T>
struct hdf5_extendable_dataset {
    _hdf5_extendable_dataset base;

    //
    // The 'bitshuffle' argument has the following meaning:
    //   0 = no compression
    //   1 = try to compress, but if plugin fails then just write uncompressed data instead
    //   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
    //   3 = compression mandatory
    //
    // List of all filter_ids: https://www.hdfgroup.org/services/contributions.html
    // Note that the compile-time constant 'bitshuffle_id' (=32008) is defined above.
    //
    hdf5_extendable_dataset(const hdf5_group &g, const std::string &dataset_name, const std::vector<hsize_t> &chunk_shape, int axis, int bitshuffle=0) :
	base(g, dataset_name, chunk_shape, axis, hdf5_type<T>(), bitshuffle)
    { }

    void write(const T *data, const std::vector<hsize_t> &shape)
    {
	base.write(data, shape);
    }

};

}  // namespace ch_frb_io

#endif // _CH_FRB_IO_INTERNALS_HPP
