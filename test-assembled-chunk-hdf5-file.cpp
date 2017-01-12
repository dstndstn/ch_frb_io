#include <cassert>
#include <iostream>
#include "ch_frb_io.hpp"

#include "assembled_chunk_msgpack.hpp"


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

    uint64_t beam_id = 42;
    int nupfreq = 16;
    int nt_per_packet = 16;
    int fpga_counts_per_sample = 400;
    uint64_t ichunk = 1000;

    std::random_device rd;
    std::mt19937 rng(rd());

    shared_ptr<assembled_chunk> chunk = assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
    chunk->randomize(rng);

    /*
     const char *filename = "test_assembled_chunk.hdf5";
     chunk->write_hdf5_file(string(filename));
     */

    chunk->msgpack_bitshuffle = true;
    string fn = "test_assembled_chunk.msgpack";
    chunk->write_msgpack_file(fn);
    cout << "Wrote to " << fn << endl;

    shared_ptr<assembled_chunk> inchunk = assembled_chunk::read_msgpack_file(fn);
    cout << "Read " << inchunk << endl;
    if (memcmp(inchunk->data, chunk->data, chunk->ndata)) {
        cout << "MISMATCH in data" << endl;
    }

    // Now fill with constant values and see that it compresses
    memset(chunk->data, 42, chunk->ndata);

    fn = "test_assembled_chunk_2.msgpack";
    chunk->write_msgpack_file(fn);
    cout << "Wrote to " << fn << endl;

    shared_ptr<assembled_chunk> inchunk2 = assembled_chunk::read_msgpack_file(fn);
    cout << "Read " << inchunk2 << endl;
    if (memcmp(inchunk2->data, chunk->data, chunk->ndata)) {
        cout << "MISMATCH in data 2" << endl;
    }



    return 0;
}
