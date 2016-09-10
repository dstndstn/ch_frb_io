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
    int nbeams = 0;
    int nupfreq = 0;
    int nfreq_coarse_per_packet = 0;
    int nt_per_packet = 0;
    int nt_tot = 0;
    int fpga_counts_per_sample = 0;
    float wt_cutoff = 0.0;

    vector<int> beam_ids;
    vector<int> freq_ids;

    unit_test_instance(std::mt19937 &rng);
};


unit_test_instance::unit_test_instance(std::mt19937 &rng)
{
    const int nfreq_coarse = ch_frb_io::constants::nfreq_coarse;

    this->nbeams = randint(rng, 1, 9);
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
	    beam_ids[i] = randint(rng, 0, constants::max_allowed_beam_id);
	} while (array_contains(&beam_ids[0], i, beam_ids[i]));
    }

    // Randomly permute frequencies, just to strengthen the test
    this->freq_ids.resize(nfreq_coarse);
    for (int i = 0; i < nfreq_coarse; i++)
	freq_ids[i] = i;
    std::shuffle(freq_ids.begin(), freq_ids.end(), rng);
}


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    auto test_instance = make_shared<unit_test_instance> (rng);
    return 0;
}
