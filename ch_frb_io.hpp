#ifndef _CH_FRB_IO_HPP
#define _CH_FRB_IO_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <string>
#include <vector>
#include <memory>
#include <random>

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

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
struct udp_packet_list;
struct udp_packet_ringbuf;
class assembled_chunk_ringbuf;


// -------------------------------------------------------------------------------------------------
//
// Compile-time constants
// Note: many of these don't really need to be fixed at compile time, and could be runtime parameters instead.


namespace constants {
    static constexpr int cache_line_size = 64;

    // Number of "coarse" (i.e. pre-upchannelized) frequency channels.
    static constexpr int nfreq_coarse_tot = 1024;

    // For an explanation of this parameter, see class intensity_network_ostream below 
    static constexpr double default_wt_cutoff = 0.3;

    // Network parameters.
    //
    // The recv_socket_timeout determines how frequently the network thread wakes up, while blocked waiting
    // for packets.  The purpose of the periodic wakeup is to check whether intensity_network_stream::end_stream()
    // has been called, and check the timeout for flushing data to assembler threads.
    //
    // The max packet size params have been chosen around ~9KB, as appropriate for 1 Gbps ethernet
    // with jumbo frames and no fragmentation.
    //
    // The stream_cancellation_latency_usec arg determines how frequently the network thread checks whether
    // intensity_network_stream::end_stream() has been called.  (We don't do this check in every iteration of
    // the packet read loop, since it requires acquiring a lock.)

    static constexpr int default_udp_port = 10252;
    static constexpr int max_input_udp_packet_size = 9000;   // largest value the input stream will accept
    static constexpr int max_output_udp_packet_size = 8910;  // largest value the output stream will produce
    static constexpr int recv_socket_timeout_usec = 10000;   // 0.01 sec
    static constexpr int stream_cancellation_latency_usec = 10000;    // 0.01 sec
    static constexpr double default_gbps = 1.0;              // default transmit rate for output stream

#ifdef __APPLE__
    // osx seems to have very small limits on socket buffer size
    static constexpr int recv_socket_bufsize = 4 * 1024 * 1024;
    static constexpr int send_socket_bufsize = 4 * 1024 * 1024;
#else
    static constexpr int recv_socket_bufsize = 128 * 1024 * 1024;
    static constexpr int send_socket_bufsize = 128 * 1024 * 1024;
#endif

    // This applies to the ring buffer between the network _output_ thread, and callers of
    // intensity_network_ostream::send_chunk().
    static constexpr int output_ringbuf_capacity = 8;

    // Ring buffer between network _input_ thread and assembler thread.
    static constexpr int unassembled_ringbuf_capacity = 16;
    static constexpr int max_unassembled_packets_per_list = 16384;
    static constexpr int max_unassembled_nbytes_per_list = 8 * 1024 * 1024;
    static constexpr int unassembled_ringbuf_timeout_usec = 250000;   // 0.25 sec

    // Ring buffers between assembler thread and processing threads.
    static constexpr int assembled_ringbuf_capacity = 8;
    static constexpr int nt_per_assembled_chunk = 1024;

    // These parameters don't really affect anything but appear in range-checking asserts.
    static constexpr int max_allowed_beam_id = 65535;
    static constexpr int max_allowed_nupfreq = 64;
    static constexpr int max_allowed_nt_per_packet = 1024;
    static constexpr int max_allowed_fpga_counts_per_sample = 1024;
};


// -------------------------------------------------------------------------------------------------
//
// intensity_network_stream: stream object which receives a packet stream from the network.
//
// (Writing a packet stream is handled by a separate class intensity_network_ostream, see below.)
//
// When the intensity_network_stream object is constructed, threads are automatically spawned to read and process
// packets.  The stream presents the incoming data to the "outside world" as a per-beam sequence 
// of regular arrays which are obtained by calling get_assembled_chunk().


class intensity_network_stream : noncopyable {
public:
    
    // The 'struct initializer' is used to construct the stream object.  A few notes:
    //
    //   - Beams are identified in the packet profile by an opaque uint16_t.  The 'beam_ids' argument
    //     should be the list of beam_ids which we expect to receive.
    //
    //   - The throw_exception_on_buffer_drop, throw_exception_on_assembler_miss flags are useful 
    //     for a unit test, but probably don't make sense otherwise.
    //
    //   - Our network protocol doesn't define any way of indicating end-of-stream.  This makes sense for the
    //     realtime search, but for testing we'd like to have a way of shutting down gracefully.  If the
    //     accept_end_of_stream_packets flag is set, then a special packet with nbeams=nupfreq=nt=0 is
    //     interpreted as an "end of stream" flag, and triggers shutdown of the network_stream.

    struct initializer {
	std::vector<int> beam_ids;
	int udp_port = constants::default_udp_port;
	bool mandate_reference_kernels = false;
	bool mandate_fast_kernels = false;
	bool emit_warning_on_buffer_drop = true;
	bool throw_exception_on_beam_id_mismatch = true;
	bool throw_exception_on_first_packet_mismatch = true;
	bool throw_exception_on_buffer_drop = false;
	bool throw_exception_on_assembler_miss = false;
	bool accept_end_of_stream_packets = true;
    };

    // Event counts are kept in an array of the form int64_t[event_type::num_types].
    // Currently we don't do anything with the event counts besides print them at the end,
    // but they're intended to be a hook for real-time monitoring via RPC's.

    enum event_type {
	byte_received = 0,
	packet_received = 1,
	packet_good = 2,
	packet_bad = 3,
	packet_dropped = 4,            // network thread will drop packets if assembler thread runs slow and ring buffer overfills
	packet_end_of_stream = 5,
	beam_id_mismatch = 6,          // beam id in packet doesn't match any element of initializer::beam_ids
	first_packet_mismatch = 7,     // stream params (nupfreq, nt_per_packet, fpga_counts_per_sample) don't match first packet received
	assembler_hit = 8,
	assembler_miss = 9,
	assembled_chunk_dropped = 10,  // assembler thread will drop assembled_chunks if processing thread runs slow
	assembled_chunk_queued = 11,
	num_types = 12                 // must be last
    };
    
    // It's convenient to initialize intensity_network_streams using a static factory function make(),
    // rather than having a public constructor.  Note that make() spawns network and assembler threads,
    // but won't listen for packets until start_stream() is called.  Between calling make() and start_stream(),
    // you'll want to spawn "consumer" threads which call get_assembled_chunk().

    static std::shared_ptr<intensity_network_stream> make(const initializer &ini_params);

    // High level control.
    void start_stream();         // tells network thread to start listening for packets
    void end_stream();           // requests stream exit (but stream will stop after a few timeouts, not immediately)
    void join_threads();         // should only be called once, does not request stream exit, blocks until network and assembler threads exit

    // This is the main routine called by the processing threads, to read data from one beam
    // corresponding to ini_params.beam_ids[assembler_ix].  (Note that the assembler_index
    // satisifes 0 <= assembler_ix < ini_params.beam_ids.size(), and is not a beam_id.)

    std::shared_ptr<assembled_chunk> get_assembled_chunk(int assembler_index);

    // Will block until first packet is received, or stream ends.  Returns true if packet was received.
    bool get_first_packet_params(int &nupfreq, int &nt_per_packet, uint64_t &fpga_counts_per_sample, uint64_t &fpga_count);
    
    // Can be called at any time, from any thread.  Note that the event counts returned by get_event_counts()
    // may slightly lag the real-time event counts (this behavior derives from wanting to avoid acquiring a
    // lock in every iteration of the packet read loop).

    initializer get_initializer();
    std::vector<int64_t> get_event_counts();

    ~intensity_network_stream();

protected:
    // Constant after construction, so not protected by lock
    const initializer ini_params;
    const int nassemblers = 0;

    // This is initialized by the assembler thread before it sets 'first_packet_received' flag.
    // Therefore, other threads can access it without a lock, but should wait for this flag to be set (which does
    // require a lock).  There is a corner case where the vector is still length-zero after the flag gets set.
    // This happens if the stream was asynchronously cancelled before receiving the first packet.

    std::vector<std::unique_ptr<assembled_chunk_ringbuf> > assemblers;
    
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

    // I'm not sure how much it actually helps bottom-line performace, but it seemed like a good idea
    // to insert padding so that data accessed by different threads is in different cache lines.
    char _pad1[constants::cache_line_size];

    // Used only by the network thread (not protected by lock)
    //
    // Note on event counting implementation: on short timescales, the network and assembler 
    // threads accumulate event counts into "local" arrays 'network_thread_event_subcounts', 
    // 'assembler_thread_event_subcounts'.  On longer timescales, these local subcounts are 
    // accumulated into the global 'cumulative_event_counts'.  This two-level accumulation
    // scheme is designed to avoid excessive contention for the event_count_lock.

    int sockfd = -1;
    std::unique_ptr<udp_packet_list> incoming_packet_list;
    std::vector<int64_t> network_thread_event_subcounts;
    char _pad2[constants::cache_line_size];

    // Used only by the assembler thread
    std::vector<int64_t> assembler_thread_event_subcounts;

    pthread_t network_thread;
    pthread_t assembler_thread;
    char _pad3[constants::cache_line_size];

    // State model.  These flags are protected by the state_lock and are set in sequence.
    // Note that the 'stream_ended_requested' flag means that the stream shutdown is imminent,
    // but doesn't mean that it has actually shut down yet, it may still be reading packets.
    // So far it hasn't been necessary to include a 'stream_ended' flag in the state model.

    pthread_mutex_t state_lock;
    pthread_cond_t cond_state_changed;       // threads wait here for state to change
    bool assembler_thread_started = false;   // set by assembler thread on startup
    bool network_thread_started = false;     // set by network thread on startup
    bool stream_started = false;             // set asynchonously by calling start_stream()
    bool first_packet_received = false;      // set by network thread
    bool assemblers_initialized = false;     // set by assembler thread
    bool stream_end_requested = false;       // can be set asynchronously by calling end_stream(), or by network/assembler threads on exit
    bool join_called = false;                // set by calling join_threads()
    char _pad4[constants::cache_line_size];

    pthread_mutex_t event_lock;
    std::vector<int64_t> cumulative_event_counts;

    // The actual constructor is protected, so it can be a helper function 
    // for intensity_network_stream::make(), but can't be called otherwise.
    intensity_network_stream(const initializer &x);

    void _open_socket();
    void _add_event_counts(std::vector<int64_t> &event_subcounts);

    static void *network_pthread_main(void *);
    static void *assembler_pthread_main(void *);

    // Private methods called by the network thread.    
    void _network_thread_body();
    void _network_thread_exit();
    void _put_unassembled_packets();

    // Private methods called by the assembler thread.     
    void _assembler_thread_body();
    void _assembler_thread_exit();
};


// -------------------------------------------------------------------------------------------------
//
// assembled_chunk
//
// This is the data structure which is returned by intensity_network_stream::get_assembled_chunk().
// It repesents a regular subarray of the intensity data in a single beam.  
//
// The data in the assembled_chunk is represented as 8-bit integers, with a coarser array of additive
// and mutiplicative offsets, but can be "decoded" to a simple floating-point array by calling
// assembled_chunk::decode().


struct assembled_chunk : noncopyable {
    // Stream parameters, specified at construction.
    const int beam_id = 0;
    const int nupfreq = 0;
    const int nt_per_packet = 0;
    const int fpga_counts_per_sample = 0;

    // More parameters which are constant after construction.
    const int nt_coarse = 0;   // equal to (constants::nt_per_assembled_chunk / nt_per_packet)
    const int nscales = 0;     // equal to (constants::nfreq_coarse * nt_coarse)
    const int ndata = 0;       // equal to (constants::nfreq_coarse * nupfreq * constants::nt_per_assembled_chunk)

    // Chunks are indexed by 'ichunk', which differs by 1 in adjacent chunks.
    //
    // Reminder: there are several units of time used in different parts of the CHIME backend!
    //   1 assembled_chunk = constants::nt_per_assembled_chunk * (1 intensity sample)
    //   1 intensity sample = assembled_chunk::fpga_counts_per_sample * (1 fpga count)
    //   1 fpga count = 2.56e-6 seconds    (approx)

    uint64_t ichunk = 0;
    uint64_t isample = 0;   // always equal to ichunk * constants::nt_per_assembled_chunk

    float *scales = nullptr;   // shape (constants::nfreq_coarse, nt_coarse)
    float *offsets = nullptr;  // shape (constants::nfreq_coarse, nt_coarse)
    uint8_t *data = nullptr;   // shape (constants::nfreq_coarse, nupfreq, constants::nt_per_assembled_chunk)

    assembled_chunk(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t ichunk);
    virtual ~assembled_chunk();
    
    // These are virtual so that subclasses can be written with optimized implementations 
    // for specific parameter choices (e.g. full CHIME nt_per_packet=16)
    virtual void add_packet(const intensity_packet &p);
    virtual void decode(float *intensity, float *weights, int stride) const;

    // Static factory function which returns either the assembled_chunk base class, or the fast_assembled_chunk
    // subclass (see below), based on the packet parameters.
    static std::shared_ptr<assembled_chunk> make(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t ichunk);

    // Utility functions currently used only for testing.
    void fill_with_copy(const std::shared_ptr<assembled_chunk> &x);
    void randomize(std::mt19937 &rng);
};


// For some choices of packet parameters (the precise criterion is nt_per_packet == 16 and
// nupfreq % 2 == 0) we can speed up assembled_chunk::add_packet() and assembled_chunk::decode()
// using assembly language kernels.

struct fast_assembled_chunk : public assembled_chunk
{
    // Constructor throws exception unless nt_per_packet == 16 and (nupfreq % 2) == 0.
    fast_assembled_chunk(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t ichunk);

    // Overrides assembled_chunk::add_packet() and assembled_chunk::decode() with fast assembly language versions.
    virtual void add_packet(const intensity_packet &p) override;
    virtual void decode(float *intensity, float *weights, int stride) const override;
};



// -------------------------------------------------------------------------------------------------
//
// intensity_network_ostream: this class is used to packetize intensity data and send
// it over the network.
//
// The stream object presents the "oustide world" with an object which is fed regular
// chunks of data via a member function send_chunk().  Under the hood, send_chunk() is
// encoding each chunk as multiple UDP packets, which are handed off to a network
// thread running in the background.


class intensity_network_ostream : noncopyable {
public:

    // The 'struct initializer' is used to construct the ostream object.  A few notes:
    //
    //   - It's convenient for the input to intensity_network_ostream to be a pair of
    //     floating-point arrays (intensity, weights).  Since the packet protocol allows
    //     a boolean mask but not floating-point weights, we simply mask all samples
    //     whose weight is below some cutoff.  This is the 'wt_cutoff' parameter.
    //
    //   - Each UDP packet is a 4D logical array of shape 
    //        (nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet). 
    //     Each chunk passed to send_chunk() is a 4D logical array of shape
    //        (nbeams, coarse_freq_ids.size(), nupfreq, nt_per_chunk)
    //     and will generally correspond to multiple packets.
    //
    //   - The target_gbps arg enables throttling of the output to some target bandwidth in Gbps.

    struct initializer {
	std::string dstname;
	std::vector<int> beam_ids;
	std::vector<int> coarse_freq_ids;   // will usually be [ 0, 1, ..., 1023 ], but reordering is allowed

	int nupfreq = 0;
	int nt_per_chunk = 0;
	int nfreq_coarse_per_packet = 0;
	int nt_per_packet = 0;
	int fpga_counts_per_sample = 0;
	float wt_cutoff = constants::default_wt_cutoff;
	double target_gbps = constants::default_gbps;   // if 0.0, then data will be written as quickly as possible!

	bool is_blocking = true;
	bool emit_warning_on_buffer_drop = true;
	bool throw_exception_on_buffer_drop = false;
    };

    const initializer ini_params;
    
    // Some of these fields are redundant with fields in the ini_params, but the redundancy is convenient.
    const int nbeams;
    const int nfreq_coarse_per_packet;
    const int nfreq_coarse_per_chunk;
    const int nupfreq;
    const int nt_per_packet;
    const int nt_per_chunk;
    const int nbytes_per_packet;
    const int npackets_per_chunk;
    const int nbytes_per_chunk;
    const int elts_per_chunk;
    const uint64_t fpga_counts_per_sample;
    const uint64_t fpga_counts_per_packet;
    const uint64_t fpga_counts_per_chunk;
    const double target_gbps;

    
    // It's convenient to initialize intensity_network_ostreams using a static factory function make(),
    // rather than having a public constructor.  Note that make() spawns a network thread which runs
    // in the background.

    static std::shared_ptr<intensity_network_ostream> make(const initializer &ini_params);

    // Data is added to the network stream by calling send_chunk().  This routine packetizes
    // the data and puts the packets in a thread-safe ring buffer, for the network thread to send.
    //
    // The 'intensity' and 'weights' arrays have logical shape
    //   (nbeams, nfreq_coarse_per_chunk, nupfreq, nt_per_chunk).
    //
    // The 'stride' arg is the memory offset between time series in the arrays.
    // To write this out explicitly, the intensity sample with logical indices (b,fcoarse,u,t) has memory location
    //   intensity + ((b + fcoarse*nupfreq) * nupfreq + u) * stride + t

    void send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count);
    
    // Called when there is no more data; sends end-of-stream packets.
    void end_stream(bool join_network_thread);

    ~intensity_network_ostream();

    // This is a helper function called by send_chunk(), but we make it public so that the unit tests can call it.
    void _encode_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count, const std::unique_ptr<udp_packet_list> &out);
    
protected:
    std::vector<uint16_t> beam_ids_16bit;
    std::vector<uint16_t> coarse_freq_ids_16bit;

    int sockfd = -1;
    std::string hostname;
    uint16_t udp_port = constants::default_udp_port;
    
    // Currently, this data is not protected by a lock, since it's only accessed by the network thread.
    // If we want to make this data accessible by other threads, this needs to be changed.
    int64_t curr_timestamp = 0;    // microseconds between first packet and most recent packet
    int64_t npackets_sent = 0;
    int64_t nbytes_sent = 0;

    // State model.
    pthread_t network_thread;
    pthread_mutex_t state_lock;
    pthread_cond_t cond_state_changed;
    bool network_thread_started = false;
    bool network_thread_joined = false;

    // Ring buffer used to exchange packets with network thread
    std::unique_ptr<udp_packet_ringbuf> ringbuf;
    
    // Temp buffer for packet encoding
    std::unique_ptr<udp_packet_list> tmp_packet_list;

    // The actual constructor is protected, so it can be a helper function 
    // for intensity_network_ostream::make(), but can't be called otherwise.
    intensity_network_ostream(const initializer &ini_params);

    static void *network_pthread_main(void *opaque_args);

    void _network_thread_body();

    void _open_socket();
    void _announce_end_of_stream();
};


// -------------------------------------------------------------------------------------------------
//
// Miscellaneous


// Utility routine: converts a string to type T (only a few T's are defined; see lexical_cast.cpp)
// Returns true on success, false on failure
template<typename T> extern bool lexical_cast(const std::string &x, T &ret);

// Also defined in lexical_cast.cpp (for the same values of T)
template<typename T> extern const char *typestr();

// Version of lexical_cast() which throws exception on failure.
template<typename T> inline T lexical_cast(const std::string &x)
{
    T ret;
    if (lexical_cast(x, ret))
	return ret;
    throw std::runtime_error("couldn't convert string '" + x + "' to " + typestr<T>());
}

// Unit tests
extern void test_lexical_cast();
extern void test_packet_offsets(std::mt19937 &rng);
extern void test_avx2_kernels(std::mt19937 &rng);
extern void peek_at_unpack_kernel();


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_HPP
