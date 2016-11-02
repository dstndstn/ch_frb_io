#include <cassert>
#include <iostream>
#include "ch_frb_io.hpp"

using namespace std;
using namespace ch_frb_io;

int main(int argc, char **argv)
{
    /*
     cout << "Warning: this will create a ~100 MB file in the current directory!\n"
     << "If this is OK press return.  If not, press control-C!\n"
     << "I AM WAITING, HUMAN: "
     << flush;
    string dummy;
    getline(cin, dummy);
     */

    const char *filename = "test_assembled_chunk.hdf5";

    uint64_t beam_id = 42;
    int nupfreq = 16;
    int nt_per_packet = 16;
    int fpga_counts_per_sample = 400;
    uint64_t ichunk = 1000;

    std::random_device rd;
    std::mt19937 rng(rd());

    shared_ptr<assembled_chunk> chunk = assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
    chunk->randomize(rng);

    chunk->write_hdf5_file(string(filename));

    
    return 0;
}
