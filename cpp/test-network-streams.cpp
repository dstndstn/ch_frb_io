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


// -------------------------------------------------------------------------------------------------


struct unit_test_instance {
    static constexpr int maxbeams = 8;

    int nbeams = 0;
    int nupfreq = 0;
    int nfreq_coarse_per_packet = 0;
    int nt_per_packet = 0;
    int nt_tot = 0;
    int fpga_counts_per_sample = 0;
    float wt_cutoff = 0.0;

    vector<int> beam_ids;
    vector<int> freq_ids;

    // not protected by lock
    pthread_t consumer_threads[maxbeams];
    
    shared_ptr<ch_frb_io::intensity_network_stream> istream;
    shared_ptr<ch_frb_io::intensity_network_ostream> ostream;    

    pthread_mutex_t lock;

    unit_test_instance(std::mt19937 &rng);
    ~unit_test_instance();
};


unit_test_instance::unit_test_instance(std::mt19937 &rng)
{
    const int nfreq_coarse = ch_frb_io::constants::nfreq_coarse;

    this->nbeams = randint(rng, 1, maxbeams+1);
    this->nupfreq = randint(rng, 1, 17);
    this->nfreq_coarse_per_packet = 1 << randint(rng,0,5);

    int n3 = nbeams * nfreq_coarse_per_packet * nupfreq;
    int nt_max = min(512, (8192+n3-1)/n3);

    this->nt_per_packet = randint(rng, nt_max/2, nt_max+1);
    this->nt_tot = nt_per_packet * randint(rng, 5000, 10000);

    this->fpga_counts_per_sample = randint(rng, 1, 1025);
    this->wt_cutoff = uniform_rand(rng, 0.3, 0.7);

    // Clunky way of generating random beam_ids
    this->beam_ids.resize(nbeams);
    for (int i = 0; i < nbeams; i++) {
	do {
	    beam_ids[i] = randint(rng, 0, ch_frb_io::constants::max_allowed_beam_id);
	} while (array_contains(&beam_ids[0], i, beam_ids[i]));
    }

    // Randomly permute frequencies, just to strengthen the test
    this->freq_ids.resize(nfreq_coarse);
    for (int i = 0; i < nfreq_coarse; i++)
	freq_ids[i] = i;
    std::shuffle(freq_ids.begin(), freq_ids.end(), rng);

    xpthread_mutex_init(&this->lock);
}


unit_test_instance::~unit_test_instance()
{
    pthread_mutex_destroy(&lock);
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


static void spawn_all_receiving_threads(const shared_ptr<unit_test_instance> &tp, std::mt19937 &rng)
{
    // Permute the beams relative to sender, to slightly strengthen the test
    vector<int> beam_ids = tp->beam_ids;
    std::shuffle(beam_ids.begin(), beam_ids.end(), rng);

    vector<shared_ptr<ch_frb_io::intensity_beam_assembler> > assemblers;

    for (int ibeam = 0; ibeam < tp->nbeams; ibeam++) {
	auto assembler = ch_frb_io::intensity_beam_assembler::make(beam_ids[ibeam]);
	spawn_consumer_thread(assembler, tp, ibeam);
	assemblers.push_back(assembler);
    }

    tp->istream = intensity_network_stream::make(assemblers, ch_frb_io::constants::default_udp_port);
    // tp->istream->start_stream();
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    for (int iouter = 0; iouter < 100; iouter++) {
	auto tp = make_shared<unit_test_instance> (rng);
	spawn_all_receiving_threads(tp, rng);

	tp->istream->end_stream(true);  // join_threads=true

	for (int ibeam = 0; ibeam < tp->nbeams; ibeam++) {
	    int err = pthread_join(tp->consumer_threads[ibeam], NULL);
	    if (err)
		throw runtime_error("pthread_join() failed");
	}
    }

    return 0;
}
