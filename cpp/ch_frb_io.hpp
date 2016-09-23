#ifndef _CH_FRB_IO_HPP
#define _CH_FRB_IO_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <string>
#include <vector>
#include <memory>
#include <random>
#include <hdf5.h>

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T> struct hdf5_extendable_dataset;

struct noncopyable
{
    noncopyable() { }
    noncopyable(const noncopyable &) = delete;
    noncopyable& operator=(const noncopyable &) = delete;
};

// Defined later in this file
struct assembled_chunk;

// Defined in ch_frb_io_internals.hpp
struct intensity_packet;
struct udp_packet_ringbuf;
class assembled_chunk_ringbuf;


// -------------------------------------------------------------------------------------------------
//
// Compile-time constants


namespace constants {
    // Number of "coarse" (i.e. pre-upchannelized) frequency channels.
    static constexpr int nfreq_coarse = 1024;

    //
    // Network parameters.
    //
    // The recv_socket_timeout is so we can periodically check for RPC's, flush data to assembler
    // threads, and notice if intensity_network_ostream::end_stream() is called.
    //
    // max_output_udp_packet_size: largest packet the output stream will produce
    // max_input_udp_packet_size: largest packet the input stream will accept (should be larger of the two)
    // We choose values around ~9KB, as appropriate for 1 Gbps ethernet with jumbo frames and no fragmentation.
    //
    // max_gbps_for_testing: This is a very small value (0.1 Gbps) but seems to be necessary if we want to 
    // avoid dropping packets in the unit test.  This means that the unit tests take about an hour to run,
    // which isn't really a problem, but is it indicative of deeper performance problems?  It would be nice
    // to understand where the bottleneck is.  (FIXME?)
    //
    static constexpr int default_udp_port = 10252;
    static constexpr int max_input_udp_packet_size = 9000;
    static constexpr int max_output_udp_packet_size = 8910;
    static constexpr int recv_socket_timeout_usec = 10000;  // 0.01 sec
    static constexpr int recv_socket_bufsize = (1 << 22);   // 4 MB
    static constexpr int send_socket_bufsize = (1 << 22);   // 4 MB
    static constexpr double max_gbps_for_testing = 0.1;
    static constexpr int stream_cancellation_latency_usec = 10000;    // 0.01 sec

    // Parameters of ring buffer between output stream object and network output thread
    static constexpr int output_ringbuf_capacity = 16;

    // Parameters of ring buffer between network thread and assembler thread
    // FIXME this seems excessive.
    static constexpr int unassembled_ringbuf_capacity = 64;
    static constexpr int max_unassembled_packets_per_list = 16384;
    static constexpr int max_unassembled_nbytes_per_list = 8 * 1024 * 1024;
    static constexpr int unassembled_ringbuf_timeout_usec = 250000;   // 0.25 sec

    // Parameters of ring buffers between assembler thread and pipeline threads.
    static constexpr int assembled_ringbuf_capacity = 8;
    static constexpr int nt_per_assembled_chunk = 1024;
    static constexpr int nt_assembler = 2 * nt_per_assembled_chunk;

    // These parameters don't really affect anything but appear in range-checking asserts.
    static constexpr int max_allowed_beam_id = 65535;
    static constexpr int max_allowed_nupfreq = 64;
    static constexpr int max_allowed_fpga_counts_per_sample = 1024;
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
// Helper classes for network code


struct udp_packet_list {
    // Initialized at construction.
    int max_npackets = 0;
    int max_nbytes = 0;
    std::unique_ptr<uint8_t[]> buf;   // points to an array of length (max_nbytes + max_packet_size).
    std::unique_ptr<int[]> off_buf;   // points to an array of length (max_npackets + 1).

    // Current state of buffer.
    int curr_npackets = 0;
    int curr_nbytes = 0;   // total size of all packets
    bool is_full = true;

    // Bare pointers.
    uint8_t *data_start = nullptr;    // points to &buf[0]
    uint8_t *data_end = nullptr;      // points to &buf[nbyes]
    int *packet_offsets = nullptr;    // points to &off_buf[0].  Note that packet_offsets[npackets] is always equal to 'nbytes'.

    udp_packet_list() { }   // default constructor makes a dummy udp_packet_list with max_npackets=max_nbytes=0.
    udp_packet_list(int max_npackets, int max_nbytes);

    // To add a packet, first add data by hand at 'data_end', then call add_packet().
    void add_packet(int packet_nbytes);

    // Doesn't deallocate buffers or change the max_* fields, but sets the current packet list to zero.
    void reset();
};


// -------------------------------------------------------------------------------------------------
//
// Network ostream


//
// intensity_network_ostream: this class is used to packetize intensity data and send
// it over the network.
//  
class intensity_network_ostream : noncopyable {
public:
    const std::string dstname;
    const int nbeams = 0;
    const int nfreq_coarse_per_chunk = 0;
    const int nfreq_coarse_per_packet = 0;
    const int nupfreq = 0;
    const int nt_per_chunk = 0;
    const int nt_per_packet = 0;
    const int fpga_counts_per_sample = 0;
    const int nbytes_per_packet = 0;
    const int npackets_per_chunk = 0;
    const int nbytes_per_chunk = 0;
    const float wt_cutoff = 0.0;
    const double target_gbps = 0.0;

    const std::vector<int> beam_ids;         // length nbeaams
    const std::vector<int> coarse_freq_ids;  // length nfreq_coarse_per_chunk

    //
    // This factory function is the "de facto constructor", used to create a new intensity_network_ostream.
    // When the intensity_network_ostream is created, a network thread is automatically spawned, which runs
    // in the background.  Data is added to the network stream by calling send_chunk().  This routine packetizes
    // the data and puts the packets in a thread-safe ring buffer, for the network thread to send.
    //
    // If the 'target_gbps' argument is nonzero, then output will be "throttled" to the target bandwidth, specified
    // in Gbps.  If target_gbps=0, then packets will be sent as quickly as possible.
    //
    static auto make(const std::string &dstname, const std::vector<int> &beam_ids,
		     const std::vector<int> &coarse_freq_ids, int nupfreq, int nt_per_chunk,
		     int nfreq_coarse_per_packet, int nt_per_packet, int fpga_counts_per_sample, 
		     float wt_cutoff, double target_gbps) 
	-> std::shared_ptr<intensity_network_ostream>;

    //
    // Called from 'external' context (i.e. same context that called make())
    // 
    // The 'intensity' and 'weights' arrays have logical shape
    //   (nbeam, nfreq_per_chunk, nupfreq, nt_per_chunk).
    //
    // The 'stride' arg is the memory offset between time series whose (beam, freq_coarse, upfreq) are consecutive.
    //
    void send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count, bool is_blocking=true);
    
    // Can be called from either external context or network thread
    void end_stream(bool join_network_thread);

    ~intensity_network_ostream();
    

protected:
    int sockfd = -1;

    std::string hostname;
    uint16_t udp_port = constants::default_udp_port;

    std::vector<uint16_t> beam_ids_16bit;
    std::vector<uint16_t> coarse_freq_ids_16bit;
    
    // Currently, this data is not protected by a lock, since it's only accessed by the network thread.
    // If we want to make this data accessible by other threads, this needs to be changed.
    int64_t curr_timestamp = 0;    // microseconds between first packet and most recent packet
    int64_t npackets_sent = 0;
    int64_t nbytes_sent = 0;

    pthread_t network_thread;
    pthread_mutex_t state_lock;
    pthread_cond_t cond_state_changed;
    bool network_thread_started = false;
    bool network_thread_joined = false;

    // ring buffer used to exchange packets with network thread
    std::unique_ptr<udp_packet_ringbuf> ringbuf;
    
    // Buffers for packet encoding
    udp_packet_list tmp_packet_list;

    // Real constructor is protected
    intensity_network_ostream(const std::string &dstname, const std::vector<int> &beam_ids, 
			      const std::vector<int> &coarse_freq_ids, int nupfreq, int nt_per_chunk,
			      int nfreq_coarse_per_packet, int nt_per_packet, int fpga_counts_per_sample, 
			      float wt_cutoff, double target_gbps);

    static void *network_pthread_main(void *opaque_args);

    void _network_thread_start();
    void _network_thread_body();

    void _open_socket();
    void _announce_end_of_stream();
};


// -------------------------------------------------------------------------------------------------
//
// The packet input stream case is more complicated!


class intensity_network_stream : noncopyable {
public:
    // Note: If 'drops_allowed' is false, the process will crash if any ring buffer overfills.
    // This sometimes makes sense during testing but probably not otherwise.
    struct initializer {
	std::vector<int> beam_ids;
	int udp_port = constants::default_udp_port;
	bool mandate_reference_kernels = false;
	bool mandate_fast_kernels = false;
	bool drops_allowed = true;
    };
    
    // It's convenient to initialize intensity_network_streams using a static factory function make(),
    // rather than having a public constructor.  Note that make() spawns a network and assembler thread,
    // but won't listen for packets until start_stream() is called.  Between calling make() and start_stream(),
    // you'll want to spawn "consumer" threads which query the assemblers for intensity and weights arrays.

    static std::shared_ptr<intensity_network_stream> make(const initializer &x);

    // High level control.
    void start_stream();         // tells network thread to start listening for packets
    void end_stream();           // asynchronously stops stream
    void join_threads();         // should only be called once, does not asynchronously stop stream.

    // FIXME Not sure if we really need this.
    // bool wait_for_first_packet(int &nupfreq, int &nt_per_packet, int &fpga_counts_per_sample);

    std::shared_ptr<assembled_chunk> get_assembled_chunk(int assembler_index);
    
    // Can be called at any time, from any thread.
    initializer get_initializer() const;
    std::vector<int64_t> get_event_counts() const;

    enum event_type {
	packet_received = 0,
	packet_good = 1,
	packet_bad = 2,
	packet_dropped = 3,
	packet_end_of_stream = 4,
	beam_id_mismatch = 5,
	first_packet_mismatch = 6,
	assembler_hits = 7,
	assembler_misses = 8,
	assembled_chunk_dropped = 9,
	assembled_chunk_processed = 10,
	num_types = 11
    };

    ~intensity_network_stream();

private:
    // FIXME do we need this?
    friend class assembled_chunk_ringbuf;

    // Constant after construction, so not protected by lock
    const initializer _initializer;
    const int nassemblers = 0;

    // This is initialized by the assembler thread before it sets 'first_packet_received' flag.
    // Therefore, other threads can access it without a lock, but should wait for this flag to be set (which does
    // require a lock).  There is a corner case where the vector is still length-zero after the flag gets set.
    // This happens if the stream was asynchronously cancelled before receiving the first packet.

    std::vector<std::shared_ptr<assembled_chunk_ringbuf> > assemblers;
    
    // These fields are initialized from the first packet received ("fp_" stands for "first packet").
    // They are initialized by the assembler thread, which then advances the state model to "first_packet_received".
    // They are not protected by a lock!  This is OK as long as non-assembler threads access them read-only, and
    // only after checking the first_packet_received flag (which does require a lock).

    uint16_t fp_nupfreq = 0;
    uint16_t fp_nt_per_packet = 0;
    uint16_t fp_fpga_counts_per_sample = 0;
    uint64_t fp_fpga_count = 0;

    // Used to exchange data between the network and assembler threads
    std::unique_ptr<udp_packet_ringbuf> unassembled_ringbuf;

    // Used only by the network thread (not protected by lock)
    int sockfd = -1;
    udp_packet_list incoming_packet_list;

    pthread_t network_thread;
    pthread_t assembler_thread;

    mutable pthread_mutex_t lock;

    // State model.  [ XXX comment on meaning of stream_ended ]
    bool assembler_thread_started = false;
    bool network_thread_started = false;
    bool stream_started = false;
    bool first_packet_received = false;
    bool stream_ended = false;
    bool join_called = false;
    pthread_cond_t cond_state_changed;

    std::vector<int64_t> event_counts;
    std::vector<int64_t> network_thread_event_subcounts;
    std::vector<int64_t> assembler_thread_event_subcounts;

    // The actual constructor is private, so it can be a helper function 
    // for intensity_network_stream::make(), but can't be called otherwise.
    intensity_network_stream(const initializer &x);

    void _open_socket();
    void _add_event_counts(std::vector<int64_t> &event_subcounts);

    static void *network_pthread_main(void *);
    static void *assembler_pthread_main(void *);

    // Private methods called by the network thread.    
    void _network_thread_start();
    void _network_thread_body();
    void _network_thread_exit();
    void _put_unassembled_packets();

    // Private methods called by the assembler thread.     
    void _assembler_thread_start();
    void _assembler_thread_body();
    void _assembler_thread_exit();
};


struct assembled_chunk : noncopyable {
    // Stream parameters
    const int beam_id = 0;
    const int nupfreq = 0;
    const int nt_per_packet = 0;
    const int fpga_counts_per_sample = 0;

    // More parameters which are constant after construction.
    const int nt_coarse = 0;   // equal to (constants::nt_per_assembled_chunk / nt_per_packet)
    const int nscales = 0;     // equal to (constants::nfreq_coarse * nt_coarse)
    const int ndata = 0;       // equal to (constants::nfreq_coarse * nupfreq * constants::nt_per_assembled_chunk)

    // Time index of first sample in chunk.
    uint64_t chunk_t0 = 0;
    uint64_t chunk_t1 = 0;

    float *scales = nullptr;   // shape (constants::nfreq_coarse, nt_coarse)
    float *offsets = nullptr;  // shape (constants::nfreq_coarse, nt_coarse)
    uint8_t *data = nullptr;   // shape (constants::nfreq_coarse, nupfreq, constants::nt_per_assembled_chunk)

    assembled_chunk(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t chunk_t0);
    virtual ~assembled_chunk();
    
    // These are virtual so that subclasses can be written with optimized implementations 
    // for specific parameter choices (e.g. full CHIME nt_per_packet=16)
    virtual void add_packet(const intensity_packet &p);
    virtual void decode(float *intensity, float *weights, int stride) const;

    // Factory function which returns either an instance of the assembled_chunk base class, or one of its subclasses.
    static std::shared_ptr<assembled_chunk> make(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t chunk_t0);

    // Utility functions currently used only for testing.
    void fill_with_copy(const std::shared_ptr<assembled_chunk> &x);
    void randomize(std::mt19937 &rng);
};


// Special case nt_per_packet=16 optimized with avx2 kernels, used in full CHIME
struct fast_assembled_chunk : public assembled_chunk
{
    fast_assembled_chunk(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t chunk_t0);

    virtual void add_packet(const intensity_packet &p) override;
    virtual void decode(float *intensity, float *weights, int stride) const override;
};


// -------------------------------------------------------------------------------------------------


// Utility routine: converts a string to type T (only a few T's are defined; see lexical_cast.cpp)
template<typename T> extern T lexical_cast(const std::string &x);

// Unit tests
extern void test_lexical_cast();
extern void test_packet_encoding();
extern void test_fast_decode_kernel(std::mt19937 &rng);
extern void peek_at_unpack_kernel();


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_HPP
