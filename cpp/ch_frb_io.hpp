#ifndef _CH_FRB_IO_HPP
#define _CH_FRB_IO_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <string>
#include <vector>
#include <memory>
#include <hdf5.h>

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T> struct hdf5_extendable_dataset;

// helper class defined in intensity_network_ostream.cpp
struct chunk_exchanger;


struct noncopyable
{
    noncopyable() { }
    noncopyable(const noncopyable &) = delete;
    noncopyable& operator=(const noncopyable &) = delete;
};


// -------------------------------------------------------------------------------------------------
//
// HDF5 file I/O
//
// Note that there are two classes here: one for reading (intensity_hdf5_file),
// and one for writing (intensity_hdf5_ofile).


struct intensity_hdf5_file : noncopyable {
    std::string filename;
    
    int nfreq;
    int npol;
    int nt_file;     // number of time samples in file (which can have gaps)
    int nt_logical;  // total number of samples in time range spanned by file, including gaps

    //
    // We currently throw an exception unless the frequencies are equally spaced and consecutive.
    // Therefore, we don't keep a full list of frequencies, just the endpoints of the range and number of samples.
    //
    // This still leaves two possibilities: the frequencies can either be ordered from lowest to highest,
    // or vice versa.  The 'frequencies_are_increasing' flag is set in the former case.  (Note that the
    // standard CHIME convention is frequencies_are_increasing=false.)
    //
    // The 'freq_lo_MHz' and 'freq_hi_MHz' fields are the lowest and highest frequencies in the entire
    // band.  Just to spell this out in detail, if frequencies_are_increasing=true, then the i-th channel
    // spans the frequency range
    //
    //   [ freq_lo + i*(freq_hi-freq_lo)/nfreq, freq_lo + (i+1)*(freq_hi-freq_lo)/nfreq ]
    //
    // and if frequences_are_increasing=false, then the i-th channel spans frequency range
    //
    //   [ freq_hi - (i+1)*(freq_hi-freq_lo)/nfreq, freq_hi - i*(freq_hi-freq_lo)/nfreq ]
    //
    // where 0 <= i <= nfreq-1.
    //
    bool frequencies_are_increasing;
    double freq_lo_MHz;
    double freq_hi_MHz;
    
    //
    // We distinguish between "file" time indices, which span the range 0 <= it < nt_file,
    // and "logical" time indices, which span the range 0 <= it < nt_logical.
    //
    // The i-th logical time sample spans the time range
    //
    //    [ time_lo + i*dt_sample, time_lo + (i+1)*dt_sample ]
    //
    double dt_sample;
    double time_lo;
    double time_hi;   // always equal to time_lo + nt_logical * dt_sample

    std::vector<double> times;                // 1D array of length nt_file
    std::vector<int> time_index_mapping;      // 1D array of length nt_file, which maps a "file" index to a "logical" index.

    // Polarization info (currently read from file but not really used)
    std::vector<std::string> polarizations;  // 1D array of length npol, each element is either "XX" or "YY"

    // 3d arrays of shape (nfreq, npol, nt_file), verbatim from the hdf5 file.
    // Rather than using them directly, you may want to use the member function get_unpolarized_intensity() below.
    std::vector<float> intensity;
    std::vector<float> weights;

    // Summary statistics
    double frac_ungapped;    // fraction of "logical" time samples which aren't in time gaps
    double frac_unmasked;    // fraction of _ungapped_ data with large weight

    // Construct from file.  If 'noisy' is true, then a one-line message will be printed when the file is read.
    explicit intensity_hdf5_file(const std::string &filename, bool noisy=true);

    //
    // Extracts a 2D array containing total intensity in time range [out_t0, out_t0+out_nt),
    // summing over polarizations.  The 'out_t0' and 'out_nt' are "logical" time indices, not
    // "file" time indices.  If this range of logical time indices contains gaps, the corresponding
    // entries of the 'out_int' and 'out_wt' arrays will be filled with zeros.
    // 
    // The 'out_int' and 'out_wt' arrays have shape (nfreq, out_nt).
    //
    // The 'out_stride' arg can be negative, if reversing the channel ordering is desired.
    // If out_stride is zero, it defaults to out_nt.
    //
    void get_unpolarized_intensity(float *out_int, float *out_wt, int out_t0, int out_nt, int out_stride=0) const;

    void run_unit_tests() const;
};


struct intensity_hdf5_ofile {
    std::string filename;
    double dt_sample;
    int nfreq;
    int npol;

    ssize_t curr_nt;       // current size of file (in time samples, not including gaps)
    double curr_time;      // time in seconds relative to arbitrary origin
    ssize_t curr_ipos;     // keeps track of gaps
    ssize_t initial_ipos;

    // used internally to print summary info when file is written
    double wsum;
    double wmax;

    std::unique_ptr<hdf5_extendable_dataset<double> > time_dataset;
    std::unique_ptr<hdf5_extendable_dataset<float> > intensity_dataset;
    std::unique_ptr<hdf5_extendable_dataset<float> > weights_dataset;

    //
    // The 'pol' argument is typically { "XX", "YY" }.  Note that pol.size() determines intensity_hdf5_file::npol,
    // which in turn determines the expected shape of the 'intensity' and 'weights' arguments passed to append_chunk().
    //
    // The 'freq0_MHz' and 'freq1_MHz' args should be the edges of the band, ordered the same way as the channel
    // indices.  For example in CHIME data, frequencies are ordered from highest to lowest, so we take freq0_MHz=800
    // and freq1_MHz=400.
    //
    // The optional 'ipos0' and 'time0' args are:
    //   ipos0 = index of first sample in file (in downsampled units, i.e. one sample is ~1.3 msec, not ~2.5 usec)
    //   time0 = arrival time of first sample in file (in seconds).
    //
    // The meaning of the 'bitshuffle' arg is:
    //   0 = no compression
    //   1 = try to compress, but if plugin fails then just write uncompressed data instead
    //   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
    //   3 = compression mandatory
    //
    // The default nt_chunk=128 comes from ch_vdif_assembler chunk size, assuming downsampling by factor 512.
    //
    intensity_hdf5_ofile(const std::string &filename, int nfreq, const std::vector<std::string> &pol,
			 double freq0_MHz, double freq1_MHz, double dt_sample, ssize_t ipos0=0,
			 double time0=0.0, int bitshuffle=2, int nt_chunk=128);

    // Note that there is no close() member function.  The file is flushed to disk and closed when the
    // intensity_hdf5_ofile destructor is called.
    ~intensity_hdf5_ofile();
    
    //
    // Append a chunk of data, of length 'nt_chunk.
    // The 'intensity' and 'weight' arrays have shape (nfreq, npol, nt_chunk).
    // The mandatory 'chunk_ipos' arg is the index of the first sample in the chunk (in downsampled units).
    // The optional 'chunk_t0' arg is the arrival time of the first sample in the chunk (if omitted, will be inferred from 'chunk_ipos').
    //
    void append_chunk(ssize_t nt_chunk, float *intensity, float *weights, ssize_t chunk_ipos, double chunk_t0);
    void append_chunk(ssize_t nt_chunk, float *intensity, float *weights, ssize_t chunk_ipos);
};


// -------------------------------------------------------------------------------------------------
//
// Network ostream


//
// The ostream writes data in "chunks", which are packetized into one or more packets.
//
// The constructor spawns a network thread.
//
// FIXME implement throughput target, reordering...
//  
struct intensity_network_ostream : noncopyable {
    const int nbeam;
    const int nupfreq;
    const int nfreq_per_chunk;
    const int nfreq_per_packet;
    const int nt_per_chunk;
    const int nt_per_packet;
    const int fpga_counts_per_sample;
    const float wt_cutoff;

    const std::vector<uint16_t> ibeam;
    const std::vector<uint16_t> ifreq_chunk;

    pthread_t network_thread;
    bool network_thread_valid = false;

    // shared data structure for communicating with network thread
    std::shared_ptr<chunk_exchanger> exchanger;
    
    // buffers for packet encoding
    std::vector<float> tmp_intensity_vec;
    std::vector<float> tmp_mask_vec;
    uint8_t *packet_buf = nullptr;
    int nbytes_per_packet = 0;
    
    intensity_network_ostream(const std::string &dstname, const std::vector<int> &ibeam, 
			      const std::vector<int> &ifreq_chunk, int nupfreq, int nt_per_chunk,
			      int nfreq_per_packet, int nt_per_packet, int fpga_counts_per_sample, 
			      float wt_cutoff);

    ~intensity_network_ostream();

    void send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count);
    
    void end_stream();
};


// -------------------------------------------------------------------------------------------------
//
// The packet input stream case is more complicated!


//
// Helper class which exchanges data between intensity_packet_stream and intensity_beam_assembler.
//
// Warning: this is a "dumb" struct with no C++ constructor/destructor/copy/move ops!
// Caller is responsible for managing ownership and making sure initialize() gets paired with destroy()
//
struct udp_packet_list : noncopyable {
    static constexpr int max_packets = 512;
    static constexpr int max_bytes = 1024 * 1024;
    static constexpr int max_packet_size = 16384;

    int beam_id;
    int npackets;
    int nbytes;
    bool is_full;

    uint8_t *data_start;    // capacity=max_bytes, size=nbytes
    uint8_t *data_end;      // always points to (data_start + nbytes)
    int *packet_offsets;    // capacity=(max_packets+1), size=(npackets+1), element at index 'npackets' is equal to 'nbytes'
    
    void initialize(int beam_id);
    void destroy();

    void swap(udp_packet_list &x);
    void clear();   // doesn't deallocate buffers

    // assumes caller has already appended packet data at (data_buf + nbytes), returns true if buffer is full
    void add_packet(int nbytes);
};


struct assembled_chunk : noncopyable {
    static constexpr int nt_per_chunk = 1024;

    // Stream parameters
    int beam_id;
    int nupfreq;   // upsampling factor (number of channels is 1024 * nupfreq)
    int fpga_counts_per_sample;

    // Time index of first sample in chunk.
    // The FPGA count of the first sample is chunk_t0 * fpga_counts_per_sample.
    uint64_t chunk_t0;

    // Arrays of shape (1024 * nupfreq, nt_per_chunk)
    // FIXME hardcoded 1024 here
    float *intensity = nullptr;
    float *weights = nullptr;

    assembled_chunk(int beam_id, int nupfreq, int fpga_counts_per_sample, uint64_t chunk_t0);
    ~assembled_chunk();
};


//
// An intensity_beam_assembler object is always backed by an assembler thread, running in the background.
// An "upstream" thread must supply the assembler's input, by calling intensity_beam_assembler::put_unassembled_packets().
// A "downstream" thread must fetch the assembler's output, by calling intensity_beam_assembler::get_assembled_packets().
//
// The assembler starts its shutdown process when intensity_beam_assembler::end_assembler() is called.  This call will
// probably be made by the "upstream" thread when it reaches end-of-stream, but we don't need to assume this.  After
// end_assembler() is called, any subsequent calls to put_unassembled_packets() will throw exceptions.  Calls to
// get_assembled_packets() will succeed until all packets have been extracted, then this call will start returning
// 'false' to indicate end-of-stream.
//
class intensity_beam_assembler : noncopyable {
public:
    const int beam_id;

    // The function intensity_beam_assembler::make() is the de facto intensity_beam_assembler constructor.
    // Thus intensity_beam_assemblers are always used through shared_ptrs.  The reason we do this is that
    // the assembler thread is always running in the background and needs to keep a reference via a shared_ptr.
    static std::shared_ptr<intensity_beam_assembler> make(int beam_id);
    
    // Helper function called by intensity_beam_assembler::make()
    void wait_for_assembler_thread_startup();

    // Called by "upstream" thread.  For a description of the 'packet_list' semantics, see the .cpp file.
    void start_stream(int fpga_counts_per_sample, int nupfreq);
    bool put_unassembled_packets(udp_packet_list &packet_list);
    void end_stream(bool join_thread);   // can also be called by assembler thread, if it exits unexpectedly

    //
    // Called by assembler thread.  
    //
    // The function get_unassembled_packets() returns true on success, false if end_assembler() has been called.  
    // In the latter case, the assembler should exit.  The 'packet_list' arg  should be an unused buffer, which 
    // is swapped for the buffer of new packets.
    //
    void assembler_thread_startup();
    bool get_unassembled_packets(udp_packet_list &packet_list);
    void put_assembled_chunk(const std::shared_ptr<assembled_chunk> &chunk);

    // Called by "downstream" thread
    bool wait_for_stream_params(int &fpga_counts_per_sample, int &nupfreq);
    bool get_assembled_chunk(std::shared_ptr<assembled_chunk> &chunk);

    ~intensity_beam_assembler();

private:
    static constexpr int unassembled_ringbuf_capacity = 16;
    static constexpr int assembled_ringbuf_capacity = 16;

    // The actual constructor is private, so it can be a helper function 
    // for intensity_beam_assembler::make(), but can't be called otherwise.
    intensity_beam_assembler(int beam_id);

    // All state below is protected by a single lock (FIXME could be made more granular)
    pthread_mutex_t lock;

    pthread_t assembler_thread;

    // Assembler state model
    bool assembler_thread_started = false;
    bool stream_started = false;
    bool stream_ended = false;
    bool assembler_thread_joined = false;
    pthread_cond_t cond_assembler_state_changed;

    // Stream parameters 
    int fpga_counts_per_sample = 0;
    int nupfreq = 0;

    udp_packet_list unassembled_ringbuf[unassembled_ringbuf_capacity];
    int unassembled_ringbuf_pos = 0;
    int unassembled_ringbuf_size = 0;

    std::shared_ptr<assembled_chunk> assembled_ringbuf[assembled_ringbuf_capacity];
    int assembled_ringbuf_pos = 0;
    int assembled_ringbuf_size = 0;

    pthread_cond_t cond_unassembled_packets_added;  // assembler thread waits here for packets
    pthread_cond_t cond_assembled_chunks_added;     // downstream thread waits here for assembled chunks
};


//
// An intensity_network_stream object is always backed by an assembler thread, running in the background.
//
// It should be easy to generalize to the case of multiple network interfaces.  In this case I think it would be
// easiest to have one intensity_network_stream object, backed by multiple network threads.
//
struct intensity_network_stream : noncopyable {
public:
    const std::vector<std::shared_ptr<intensity_beam_assembler> > assemblers;
    const int udp_port;

    // De facto constructor.  A thread is spawned, but it won't start reading packets until start_stream() is called.
    static std::shared_ptr<intensity_network_stream> make(const std::vector<std::shared_ptr<intensity_beam_assembler> > &assemblers, int udp_port);
    
    // High level control routines.
    //   - end_stream() can be used either to cancel a stream which never ran, or interrupt a running stream.
    //   - wait_for_first_packet_params() return false if the stream gets cancelled before receiving packets.

    void start_stream();   // tells network thread to start reading packets
    void wait_for_network_thread_startup();
    void wait_for_end_of_stream(bool join_threads);
    void end_stream(bool join_threads);

    // Called by network thread.
    // When the network thread exits (typically when end-of-stream packet is received), it calls end_stream().
    void network_thread_startup();
    bool wait_for_start_stream();     // returns false if stream has been cancelled

    ~intensity_network_stream();

private:
    // All state below is protected by this lock.
    pthread_mutex_t lock;

    pthread_t network_thread;

    bool network_thread_started = false;    // set before intensity_network_stream::make() returns
    bool stream_started = false;            // set when intensity_network_stream::start_stream() is called
    bool stream_ended = false;              // set when "end-of-stream" packet arrives, or network thread unexpectedly exits
    bool network_thread_joined = false;     // set in wait_for_end_of_stream(), but only if join_threads flag is set
    pthread_cond_t cond_state_changed;

    // The actual constructor is private, so it can be a helper function 
    // for intensity_network_stream::make(), but can't be called otherwise.
    intensity_network_stream(const std::vector<std::shared_ptr<intensity_beam_assembler> > &assemblers, int udp_port);
};


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
    
    // Write multidimensional dataset
    template<typename T> void write_dataset(const std::string &dataset_name, const T *data, const std::vector<hsize_t> &shape)
    {
	this->_write_dataset(dataset_name, hdf5_type<T>(), reinterpret_cast<const void *> (data), shape);
    }

    // This interface is intended for small string-valued datasets.
    void write_string_dataset(const std::string &dataset_name, const std::vector<std::string> &data, const std::vector<hsize_t> &shape);
    void read_string_dataset(const std::string &dataset_name, std::vector<std::string> &data, const std::vector<hsize_t> &expected_shape) const;

    // Helpers
    void _get_attribute_shape(const std::string &attr_name, hid_t attr_id, std::vector<hsize_t> &shape) const;
    void _read_attribute(const std::string &attr_name, hid_t hdf5_type, void *out, const std::vector<hsize_t> &expected_shape) const;
    void _write_attribute(const std::string &attr_name, hid_t hdf5_type, const void *data, const std::vector<hsize_t> &shape);
    void _get_dataset_shape(const std::string &dataset_name, hid_t dataset_id, std::vector<hsize_t> &shape) const;
    void _check_dataset_shape(const std::string &dataset_name, hid_t dataset_id, const std::vector<hsize_t> &expected_shape) const;
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

#endif // _CH_FRB_IO_HPP
