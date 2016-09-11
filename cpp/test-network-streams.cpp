#include <cassert>
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
    return 0.823*beam_id + 1.319*ifreq + 0.2139*it;
}


// weights are between 0.2 and 1
inline double wtval(int beam_id, int ifreq, int it)
{
    return 0.6 + 0.4 * sin(1.328*beam_id + 2.382*ifreq + 0.883*it);
}


// -------------------------------------------------------------------------------------------------


struct unit_test_instance {
    static constexpr int maxbeams = 8;

    int nbeams = 0;
    int nupfreq = 0;
    int nfreq_coarse_per_packet = 0;
    int nt_per_packet = 0;
    int nt_per_chunk = 0;
    int nt_tot = 0;
    int fpga_counts_per_sample = 0;
    uint64_t initial_fpga_count = 0;
    float wt_cutoff = 0.0;

    vector<int> recv_beam_ids;
    vector<int> send_beam_ids;
    vector<int> send_freq_ids;
    int send_stride;

    // not protected by lock
    pthread_t consumer_threads[maxbeams];    
    shared_ptr<ch_frb_io::intensity_network_stream> istream;

    pthread_mutex_t lock;

    unit_test_instance(std::mt19937 &rng);
    ~unit_test_instance();

    void show() const;
};


unit_test_instance::unit_test_instance(std::mt19937 &rng)
{
    const int nfreq_coarse = ch_frb_io::constants::nfreq_coarse;

    this->nbeams = randint(rng, 1, maxbeams+1);
    this->nupfreq = randint(rng, 1, 17);
    this->nfreq_coarse_per_packet = 1 << randint(rng,0,5);
    
    // nt_max = max possible nt_per_packet (before max packet size is exceeded)
    int header_nbytes = 24 + 2*nbeams + 2*nfreq_coarse_per_packet + 8*nbeams*nfreq_coarse_per_packet;
    int nbytes_per_nt = nbeams * nfreq_coarse_per_packet * nupfreq;
    int nt_max = (ch_frb_io::constants::max_output_udp_packet_size - header_nbytes) / nbytes_per_nt;
    nt_max = min(512, nt_max);

    assert(nt_max >= 3);
    this->nt_per_packet = randint(rng, nt_max/2+1, nt_max+1);

    // FIXME increase nt_tot
    this->nt_per_chunk = nt_per_packet * randint(rng,1,5);
    this->nt_tot = nt_per_chunk * randint(rng,1,5);

    this->fpga_counts_per_sample = randint(rng, 1, 1025);
    this->initial_fpga_count = fpga_counts_per_sample * randint(rng, 0, 4097);
    this->wt_cutoff = uniform_rand(rng, 0.3, 0.7);

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
    this->send_freq_ids.resize(nfreq_coarse);
    for (int i = 0; i < nfreq_coarse; i++)
	send_freq_ids[i] = i;
    std::shuffle(send_freq_ids.begin(), send_freq_ids.end(), rng);

    this->send_stride = randint(rng, nt_per_chunk, 2*nt_per_chunk+1);

    // Will be needed shortly for synchronization
    xpthread_mutex_init(&this->lock);
}


unit_test_instance::~unit_test_instance()
{
    pthread_mutex_destroy(&lock);
}


void unit_test_instance::show() const
{
    cout << "nbeams=" << nbeams << endl
	 << "nupfreq=" << nupfreq << endl
	 << "nfreq_coarse_per_packet=" << nfreq_coarse_per_packet << endl
	 << "nt_per_packet=" << nt_per_packet << endl
	 << "nt_per_chunk=" << nt_per_chunk << endl
	 << "nt_tot=" << nt_tot << endl;
}


// -------------------------------------------------------------------------------------------------


struct consumer_thread_context {
    shared_ptr<ch_frb_io::intensity_beam_assembler> assembler;
    shared_ptr<unit_test_instance> tp;

    pthread_mutex_t lock;
    pthread_cond_t cond_running;
    bool is_running = false;

    consumer_thread_context(const shared_ptr<ch_frb_io::intensity_beam_assembler> &assembler_, const shared_ptr<unit_test_instance> &tp_)
	: assembler(assembler_), tp(tp_)
    {
	xpthread_mutex_init(&lock);
	xpthread_cond_init(&cond_running);
    }
};


static void *consumer_thread_main(void *opaque_arg)
{
    consumer_thread_context *context = reinterpret_cast<consumer_thread_context *> (opaque_arg);
    shared_ptr<ch_frb_io::intensity_beam_assembler> assembler = context->assembler;
    shared_ptr<unit_test_instance> tp = context->tp;

    pthread_mutex_lock(&context->lock);
    context->is_running = true;
    pthread_cond_broadcast(&context->cond_running);
    pthread_mutex_unlock(&context->lock);

    // actual logic of the consumer thread goes here
    cerr << "consumer_thread: artificially exiting\n";
    return NULL;
}


static void spawn_consumer_thread(const shared_ptr<ch_frb_io::intensity_beam_assembler> &assembler, const shared_ptr<unit_test_instance> &tp, int ibeam)
{
    consumer_thread_context context(assembler, tp);

    int err = pthread_create(&tp->consumer_threads[ibeam], NULL, consumer_thread_main, &context);
    if (err)
	throw runtime_error(string("pthread_create() failed to create consumer thread: ") + strerror(errno));

    pthread_mutex_lock(&context.lock);
    while (!context.is_running)
	pthread_cond_wait(&context.cond_running, &context.lock);
    pthread_mutex_unlock(&context.lock);
}


static void spawn_all_receive_threads(const shared_ptr<unit_test_instance> &tp)
{
    vector<shared_ptr<ch_frb_io::intensity_beam_assembler> > assemblers;

    for (int ibeam = 0; ibeam < tp->nbeams; ibeam++) {
	auto assembler = ch_frb_io::intensity_beam_assembler::make(tp->recv_beam_ids[ibeam]);
	spawn_consumer_thread(assembler, tp, ibeam);
	assemblers.push_back(assembler);
    }

    tp->istream = intensity_network_stream::make(assemblers, ch_frb_io::constants::default_udp_port);
}


// -------------------------------------------------------------------------------------------------


static void send_data(const shared_ptr<unit_test_instance> &tp)
{
    const int nbeams = tp->nbeams;
    const int nupfreq = tp->nupfreq;
    const int nfreq_coarse = ch_frb_io::constants::nfreq_coarse;
    const int nt_chunk = tp->nt_per_chunk;
    const int nchunks = tp->nt_tot / tp->nt_per_chunk;
    const int stride = tp->send_stride;
    const int s2 = nupfreq * stride;
    const int s3 = nfreq_coarse * s2;

    string dstname = "127.0.0.1:" + to_string(ch_frb_io::constants::default_udp_port);

    // spawns network thread
    auto ostream = intensity_network_ostream::make(dstname, tp->send_beam_ids, tp->send_freq_ids, tp->nupfreq,
						   tp->nt_per_chunk, tp->nfreq_coarse_per_packet, 
						   tp->nt_per_packet, tp->fpga_counts_per_sample,
						   tp->wt_cutoff, 0.0);

    vector<float> intensity(nbeams * s3, 0.0);
    vector<float> weights(nbeams * s3, 0.0);

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	int it0 = ichunk * nt_chunk;   // start of chunk

	for (int ibeam = 0; ibeam < nbeams; ibeam++) {
	    int beam_id = tp->send_beam_ids[ibeam];

	    for (int ifreq_coarse = 0; ifreq_coarse < nfreq_coarse; ifreq_coarse++) {
		int coarse_freq_id = tp->send_freq_ids[ifreq_coarse];

		for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
		    int ifreq_logical = coarse_freq_id * nupfreq + iupfreq;
		    
		    int row_start = ibeam*s3 + ifreq_coarse*s2 + iupfreq*stride;
		    float *i_row = &intensity[row_start];
		    float *w_row = &weights[row_start];

		    for (int it = 0; it < nt_chunk; it++) {
			i_row[it] = intval(beam_id, ifreq_logical, it+it0);
			w_row[it] = wtval(beam_id, ifreq_logical, it+it0);
		    }
		}
	    }
	}

	uint64_t fpga_count = tp->initial_fpga_count + ichunk * nt_chunk * tp->fpga_counts_per_sample;
	ostream->send_chunk(&intensity[0], &weights[0], stride, fpga_count, true);
    }

    // joins network thread
    ostream->end_stream(true);
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    for (int iouter = 0; iouter < 100; iouter++) {
	auto tp = make_shared<unit_test_instance> (rng);
	tp->show();

	spawn_all_receive_threads(tp);
	tp->istream->start_stream();

	send_data(tp);	
	tp->istream->end_stream(true);  // join_threads=true

	for (int ibeam = 0; ibeam < tp->nbeams; ibeam++) {
	    int err = pthread_join(tp->consumer_threads[ibeam], NULL);
	    if (err)
		throw runtime_error("pthread_join() failed");
	}
    }

    return 0;
}
