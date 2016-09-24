#include <cassert>
#include <algorithm>
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


inline bool array_contains(const int *arr, int len, int x)
{
    for (int i = 0; i < len; i++)
	if (arr[i] == x)
	    return true;
    return false;
}


inline double intval(int beam_id, int ifreq, int it)
{
    return fmod(0.823*beam_id + 1.319*ifreq + 1.023*it, 1.0);
}


// weights are between 0.2 and 1
inline double wtval(int beam_id, int ifreq, int it)
{
    return 0.2 + 0.8 * fmod(1.328*beam_id + 2.382*ifreq + 0.883*it, 1.0);
}


// -------------------------------------------------------------------------------------------------


struct unit_test_instance {
    static constexpr int maxbeams = 8;

    bool use_fast_kernels = false;
    int nbeams = 0;
    int nupfreq = 0;
    int nfreq_coarse_per_packet = 0;
    int nt_per_packet = 0;
    int nt_per_chunk = 0;
    int nt_tot = 0;
    int fpga_counts_per_sample = 0;
    uint64_t initial_t0 = 0;
    float wt_cutoff = 0.0;

    vector<int> recv_beam_ids;
    vector<int> send_beam_ids;
    vector<int> send_freq_ids;

    int send_stride = 0;
    int recv_stride = 0;

    int nbytes_per_packet = 0;
    int npackets_per_chunk = 0;

    // not protected by lock
    pthread_t consumer_threads[maxbeams];    
    shared_ptr<ch_frb_io::intensity_network_stream> istream;
    shared_ptr<ch_frb_io::intensity_network_ostream> ostream;

    pthread_mutex_t tpos_lock;
    pthread_cond_t cond_tpos_changed;
    uint64_t consumer_tpos[maxbeams];

    unit_test_instance(std::mt19937 &rng, int irun, int nrun);
    ~unit_test_instance();
};


unit_test_instance::unit_test_instance(std::mt19937 &rng, int irun, int nrun)
{
    const int nfreq_coarse_tot = ch_frb_io::constants::nfreq_coarse;
    const int nt_assembler = ch_frb_io::constants::nt_assembler;

    // In alternating iterations of the test, we choose parameters so that the "fast" kernels are used.
    this->use_fast_kernels = ((irun % 2) == 0);
    this->nbeams = randint(rng, 1, maxbeams+1);
    this->nupfreq = use_fast_kernels ? (2*randint(rng,1,9)) : randint(rng,1,17);
    this->nt_per_packet = use_fast_kernels ? 16 : (1 << randint(rng,0,5));

    // Now assign nfreq_coarse_per_packet, subject to packet size constraints.
    // The constants "c0" and "c1" are defined so that the packet size is c0 + c1 * nfreq_coarse_per_packet.
    int c0 = 24 + 2*nbeams;
    int c1 = 2 + 8*nbeams + nbeams*nupfreq*nt_per_packet;

    this->nfreq_coarse_per_packet = (ch_frb_io::constants::max_output_udp_packet_size - c0) / c1;
    this->nfreq_coarse_per_packet = min(nfreq_coarse_per_packet, ch_frb_io::constants::nfreq_coarse);
    this->nfreq_coarse_per_packet = round_down_to_power_of_two(nfreq_coarse_per_packet);

    assert(nfreq_coarse_per_packet >= 4);
    this->nfreq_coarse_per_packet /= (1 << randint(rng,0,3));
    
    // Assign nt_per_chunk.  Each chunk should be no more than 512 samples.
    this->nt_per_chunk = nt_per_packet * randint(rng, 1, 512/nt_per_packet + 1);

    // Assign nt_tot.  We require <= 1024 chunks, and <= 1 GB total (summed over all beams).
    // FIXME think about increasing the 1 GB limit.  (Watch out for 32-bit overflow!)
    int packet_nbytes = packet_size(nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet);
    int chunk_nbytes = packet_nbytes * (nfreq_coarse_tot / nfreq_coarse_per_packet) * (nt_per_chunk / nt_per_packet);
    int max_nchunks = min(1024, (1<<30) / chunk_nbytes);
    this->nt_tot = nt_per_chunk * randint(rng, 1, max_nchunks+1);

    // We now require that initial_t0 is a multiple of nt_per_packet.
    this->initial_t0 = randint(rng, 0, 4097) * nt_per_packet;
    this->fpga_counts_per_sample = randint(rng, 1, 1025);
    this->wt_cutoff = uniform_rand(rng, 0.3, 0.7);

    this->send_stride = randint(rng, nt_per_chunk, 2*nt_per_chunk+1);
    this->recv_stride = randint(rng, constants::nt_per_assembled_chunk, 2 * constants::nt_per_assembled_chunk);

#if 0
    // Sometimes it's convenient to debug a specific test case...
    this->use_fast_kernels = false;
    this->nbeams = 8;
    this->nupfreq = 4;
    this->nfreq_coarse_per_packet = 8;
    this->nt_per_packet = 8;
    this->nt_per_chunk = 488;
    this->nt_tot = 6344;
    this->initial_t0 = 15840;
    this->fpga_counts_per_sample = 862;
    this->wt_cutoff = 0.565123;
    this->send_stride = 646;
    this->recv_stride = 1424;
#endif

    // Clunky way of generating random beam_ids
    this->recv_beam_ids.resize(nbeams);
    for (int i = 0; i < nbeams; i++) {
	do {
	    recv_beam_ids[i] = randint(rng, 0, ch_frb_io::constants::max_allowed_beam_id);
	} while (array_contains(&recv_beam_ids[0], i, recv_beam_ids[i]));
    }

    // To slightly strengthen the test, we permute the receiver beams relative to sender
    this->send_beam_ids = recv_beam_ids;
    std::shuffle(send_beam_ids.begin(), send_beam_ids.end(), rng);

    // Randomly permute frequencies, just to strengthen the test
    this->send_freq_ids.resize(nfreq_coarse_tot);
    for (int i = 0; i < nfreq_coarse_tot; i++)
	send_freq_ids[i] = i;
    std::shuffle(send_freq_ids.begin(), send_freq_ids.end(), rng);

    this->nbytes_per_packet = packet_size(nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet);
    this->npackets_per_chunk = (nt_per_chunk / nt_per_packet) * (nfreq_coarse_tot / nfreq_coarse_per_packet);

    xpthread_mutex_init(&this->tpos_lock);
    xpthread_cond_init(&this->cond_tpos_changed);
    
    for (int ithread = 0; ithread < nbeams; ithread++)
	this->consumer_tpos[ithread] = initial_t0;

    cout << "\nStarting test run " << irun << "/" << nrun << endl;

    cout << "    use_fast_kernels=" << use_fast_kernels << endl
	 << "    nbeams=" << nbeams << endl
	 << "    nupfreq=" << nupfreq << endl
	 << "    nfreq_coarse_per_packet=" << nfreq_coarse_per_packet << endl
	 << "    nt_per_packet=" << nt_per_packet << endl
	 << "    nt_per_chunk=" << nt_per_chunk << endl
	 << "    nt_tot=" << nt_tot << endl
	 << "    initial_t0=" << initial_t0 << endl
	 << "    fpga_counts_per_sample=" << fpga_counts_per_sample << endl
	 << "    wt_cutoff=" << wt_cutoff << endl
	 << "    send_stride=" << send_stride << endl
	 << "    recv_stride=" << recv_stride << endl
	 << "    nbytes_per_packet=" << nbytes_per_packet << endl
	 << "    npackets_per_chunk=" << npackets_per_chunk << endl;

    // Worst-case storage requirements for unassembled ringbuf.
    int wc_nchunks = min(nt_assembler/nt_per_chunk + 1, nt_tot/nt_per_chunk);
    int wc_npackets = wc_nchunks * npackets_per_chunk;
    int wc_nbytes = wc_npackets * nbytes_per_packet;
    
    // Storage actually allocated
    int npackets_alloc = ch_frb_io::constants::unassembled_ringbuf_capacity * ch_frb_io::constants::max_unassembled_packets_per_list;
    int nbytes_alloc = ch_frb_io::constants::unassembled_ringbuf_capacity * ch_frb_io::constants::max_unassembled_nbytes_per_list;
    
    if ((npackets_alloc < wc_npackets) || (nbytes_alloc < wc_nbytes)) {
	cout << "    npackets_needed=" << wc_npackets << endl
	     << "    npackets_allocated=" << npackets_alloc << endl
	     << "    nbytes_needed=" << wc_nbytes << endl
	     << "    nbytes_alloc=" << nbytes_alloc << endl
	     << "Fatal: unassembled_packet_buf is underallocated" << endl;
	
	exit(1);
    }

    cout << endl;
}


unit_test_instance::~unit_test_instance()
{
    pthread_mutex_destroy(&tpos_lock);
    pthread_cond_destroy(&cond_tpos_changed);
}


// -------------------------------------------------------------------------------------------------


struct consumer_thread_context {
    shared_ptr<unit_test_instance> tp;
    int ithread = 0;

    pthread_mutex_t lock;
    pthread_cond_t cond_running;
    bool is_running = false;

    consumer_thread_context(const shared_ptr<unit_test_instance> &tp_, int ithread_)
	: tp(tp_), ithread(ithread_)
    {
	xpthread_mutex_init(&lock);
	xpthread_cond_init(&cond_running);
    }
};


static void *consumer_thread_main(void *opaque_arg)
{
    consumer_thread_context *context = reinterpret_cast<consumer_thread_context *> (opaque_arg);

    //
    // Note: the consumer thread startup logic works like this:
    //
    //   - parent thread puts a context struct on its stack, in spawn_consumer_thread()
    //   - parent thread calls pthread_create() to spawn consumer thread
    //   - parent thread blocks waiting for consumer thread to set context->is_running
    //   - when parent thread unblocks, spawn_consumer_thread() removes and the context struct becomes invalid
    //
    // Therefore, the consumer thread is only allowed to access the context struct _before_
    // setting context->is_running to unblock the parent thread.  The first thing we do is
    // extract all members of the context struct so we don't need to access it again.
    //
    shared_ptr<unit_test_instance> tp = context->tp;
    int ithread = context->ithread;
    
    // Now we can set context->is_running and unblock the parent thread.
    pthread_mutex_lock(&context->lock);
    context->is_running = true;
    pthread_cond_broadcast(&context->cond_running);
    pthread_mutex_unlock(&context->lock);

    int nalloc = ch_frb_io::constants::nfreq_coarse * tp->nupfreq * tp->recv_stride;
    vector<float> all_intensities(nalloc, 0.0);
    vector<float> all_weights(nalloc, 0.0);

    double wt_cutoff = tp->wt_cutoff;
    int test_t0 = tp->initial_t0;
    int test_t1 = tp->initial_t0 + tp->nt_tot;
    int beam_id = tp->recv_beam_ids[ithread];

    uint64_t tpos = 0;
    bool tpos_initialized = false;

    for (;;) {
	auto chunk = tp->istream->get_assembled_chunk(ithread);

	if (!chunk)
	    break;

	chunk->decode(&all_intensities[0], &all_weights[0], tp->recv_stride);

	assert(chunk->nupfreq == tp->nupfreq);
	assert(chunk->nt_per_packet == tp->nt_per_packet);
	assert(chunk->fpga_counts_per_sample == tp->fpga_counts_per_sample);
	assert(chunk->beam_id == beam_id);

	if (tpos_initialized)
	    assert(chunk->isample == tpos);
	else
	    assert(chunk->isample <= tp->initial_t0);

	// tpos = expected isample in the next assembled_chunk.
	tpos = chunk->isample + ch_frb_io::constants::nt_per_assembled_chunk;
	tpos_initialized = true;

	pthread_mutex_lock(&tp->tpos_lock);
	tp->consumer_tpos[ithread] = tpos;
	pthread_cond_broadcast(&tp->cond_tpos_changed);
	pthread_mutex_unlock(&tp->tpos_lock);

	int chunk_t0 = chunk->isample;

	for (int ifreq = 0; ifreq < ch_frb_io::constants::nfreq_coarse * tp->nupfreq; ifreq++) {
	    const float *int_row = &all_intensities[0] + ifreq * tp->recv_stride;
	    const float *wt_row = &all_weights[0] + ifreq * tp->recv_stride;

	    for (int it = 0; it < ch_frb_io::constants::nt_per_assembled_chunk; it++) {
		// Out of range
		if ((it+chunk_t0 < test_t0) || (it+chunk_t0 >= test_t1)) {
		    assert(wt_row[it] == 0.0);
		    continue;
		}

		double ival = intval(beam_id, ifreq, it+chunk_t0);
		double wval = wtval(beam_id, ifreq, it+chunk_t0);

		// Check intensity
		if ((wt_row[it] > 0.0) && fabs(int_row[it] - ival) > 0.021) {
		    stringstream ss;
		    ss << "Test failure in weights array: beam_id=" << beam_id << ", ifreq=" << ifreq << ", it=" << (it+chunk_t0) << "\n"
		       << "   intval(...)=" << ival << ", intensity_received=" << int_row[it] << "\n";

		    cerr << ss.str();
		    exit(1);
		}

		// Check weights
		if ((wt_row[it] == 0.0) && (wval <= 1.00001 * wt_cutoff))
		    continue;
		if ((wt_row[it] == 1.0) && (wval >= 0.99999 * wt_cutoff))
		    continue;

		// If we get here, the weights check failed
		stringstream ss;
		ss << "Test failure in weights array: beam_id=" << beam_id << ", ifreq=" << ifreq << ", it=" << (it+chunk_t0) << "\n"
		   << "   wtval(...)=" << wval << ", wt_cutoff=" << wt_cutoff << ", wt_received=" << wt_row[it] << "\n";

		if ((wt_row[it] == 0.0) && (wval > wt_cutoff))
		    ss << "A probable reason for this failure is that the receive threads can't keep up with the sender.\n"
		       << "Try decreasing target_gbps (in test_network-streams.cpp:send_data()), recompiling and trying again.\n";

		cerr << ss.str();
		exit(1);
	    }
	}

	// more chunk processing will go here
    }

    assert(tpos_initialized);
    assert(tpos >= tp->initial_t0 + tp->nt_tot);

    return NULL;
}


static void spawn_consumer_thread(const shared_ptr<unit_test_instance> &tp, int ithread)
{
    consumer_thread_context context(tp, ithread);

    int err = pthread_create(&tp->consumer_threads[ithread], NULL, consumer_thread_main, &context);
    if (err)
	throw runtime_error(string("pthread_create() failed to create consumer thread: ") + strerror(errno));

    pthread_mutex_lock(&context.lock);
    while (!context.is_running)
	pthread_cond_wait(&context.cond_running, &context.lock);
    pthread_mutex_unlock(&context.lock);
}


static void spawn_all_receive_threads(const shared_ptr<unit_test_instance> &tp)
{
    ch_frb_io::intensity_network_stream::initializer initializer;
    initializer.beam_ids = tp->recv_beam_ids;
    initializer.mandate_reference_kernels = !tp->use_fast_kernels;
    initializer.mandate_fast_kernels = tp->use_fast_kernels;
    initializer.throw_exception_on_buffer_drop = true;
    initializer.throw_exception_on_assembler_miss = true;

    tp->istream = intensity_network_stream::make(initializer);
    
    for (int ithread = 0; ithread < tp->nbeams; ithread++)
	spawn_consumer_thread(tp, ithread);
}


// -------------------------------------------------------------------------------------------------


static void send_data(const shared_ptr<unit_test_instance> &tp)
{
    intensity_network_ostream::initializer ini_params;
    ini_params.dstname = "127.0.0.1";
    ini_params.beam_ids = tp->send_beam_ids;
    ini_params.coarse_freq_ids = tp->send_freq_ids;
    ini_params.nupfreq = tp->nupfreq;
    ini_params.nt_per_chunk = tp->nt_per_chunk;
    ini_params.nfreq_coarse_per_packet = tp->nfreq_coarse_per_packet;
    ini_params.nt_per_packet = tp->nt_per_packet;
    ini_params.fpga_counts_per_sample = tp->fpga_counts_per_sample;
    ini_params.wt_cutoff = tp->wt_cutoff;

    // FIXME: currently we have to run test-network-streams.cpp at very low throughput (0.1 Gbps)
    // to avoid dropping packets.  This means that the unit tests take about an hour to run,
    // which isn't really a problem, but is indicative of deeper performance problems?  It
    // would be nice to understand where the bottleneck is.

    ini_params.target_gbps = 0.1;
    
    cout << "\nNote: target_gbps = " << ini_params.target_gbps << "\n";

    // spawns network thread
    tp->ostream = intensity_network_ostream::make(ini_params);

    const int nt_assembler = ch_frb_io::constants::nt_assembler;
    const int nfreq_coarse_tot = ch_frb_io::constants::nfreq_coarse;

    const int nbeams = tp->nbeams;
    const int nupfreq = tp->nupfreq;
    const int nt_chunk = tp->nt_per_chunk;
    const int nchunks = tp->nt_tot / tp->nt_per_chunk;
    const int stride = tp->send_stride;
    const int s2 = nupfreq * stride;
    const int s3 = nfreq_coarse_tot * s2;

    vector<float> intensity(nbeams * s3, 0.0);
    vector<float> weights(nbeams * s3, 0.0);

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	int chunk_t0 = tp->initial_t0 + ichunk * nt_chunk;   // start of chunk

	for (int ibeam = 0; ibeam < nbeams; ibeam++) {
	    int beam_id = tp->send_beam_ids[ibeam];

	    for (int ifreq_coarse = 0; ifreq_coarse < nfreq_coarse_tot; ifreq_coarse++) {
		int coarse_freq_id = tp->send_freq_ids[ifreq_coarse];

		for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
		    int ifreq_logical = coarse_freq_id * nupfreq + iupfreq;
		    
		    int row_start = ibeam*s3 + ifreq_coarse*s2 + iupfreq*stride;
		    float *i_row = &intensity[row_start];
		    float *w_row = &weights[row_start];

		    for (int it = 0; it < nt_chunk; it++) {
			i_row[it] = intval(beam_id, ifreq_logical, chunk_t0 + it);
			w_row[it] = wtval(beam_id, ifreq_logical, chunk_t0 + it);
		    }
		}
	    }
	}

	// Wait for consumer threads if necessary.
	// Note that for some choices of unit_test_instance parameters, this can test the timeout logic.
	pthread_mutex_lock(&tp->tpos_lock);
	for (int i = 0; i < nbeams; i++) {
	    while (tp->consumer_tpos[i] + nt_assembler < uint64_t(chunk_t0))
		pthread_cond_wait(&tp->cond_tpos_changed, &tp->tpos_lock);
	}
	pthread_mutex_unlock(&tp->tpos_lock);

	uint64_t fpga_count = (tp->initial_t0 + ichunk * nt_chunk) * tp->fpga_counts_per_sample;
	tp->ostream->send_chunk(&intensity[0], &weights[0], stride, fpga_count);
	cout << "sent chunk " << ichunk << "/" << nchunks << endl;
    }

    // joins network thread
    tp->ostream->end_stream(true);
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    const int nrun = 100;

#if 1
    std::random_device rd;
    std::mt19937 rng(rd());
#else
    // Sometimes it's convenient to use the same seed every time, for debugging.
    std::mt19937 rng;
#endif

    cout << "Warning: this cpu-intensive, multithreaded test will take over your machine for about an hour!\n"
	 << "If this is OK press return.  If not, press control-C!\n"
	 << "I AM WAITING, HUMAN: "
	 << flush;

    string dummy;
    getline(cin, dummy);

    for (int irun = 0; irun < nrun; irun++) {
	auto tp = make_shared<unit_test_instance> (rng, irun, nrun);

	spawn_all_receive_threads(tp);
	tp->istream->start_stream();

	send_data(tp);	
	tp->istream->join_threads();

	for (int ibeam = 0; ibeam < tp->nbeams; ibeam++) {
	    int err = pthread_join(tp->consumer_threads[ibeam], NULL);
	    if (err)
		throw runtime_error("pthread_join() failed");	    
	}

	vector<int64_t> counts = tp->istream->get_event_counts();

	typedef ch_frb_io::intensity_network_stream::event_type ev_type;
	assert(counts[ev_type::packet_bad] == 0);
	assert(counts[ev_type::packet_bad] == 0);
	assert(counts[ev_type::packet_dropped] == 0);
	assert(counts[ev_type::beam_id_mismatch] == 0);
	assert(counts[ev_type::first_packet_mismatch] == 0);
	assert(counts[ev_type::assembler_miss] == 0);
	assert(counts[ev_type::assembled_chunk_dropped] == 0);

	int expected_npackets = (tp->nt_tot / tp->nt_per_chunk) * tp->npackets_per_chunk;
	assert(counts[ev_type::packet_received] - counts[ev_type::packet_end_of_stream] == expected_npackets);
    }

    cout << "\n    ****  network test passed!! ****\n\n";

    return 0;
}
