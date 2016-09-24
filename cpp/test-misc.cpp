#include "ch_frb_io_internals.hpp"

using namespace ch_frb_io;

int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    test_lexical_cast();
    test_packet_offsets(rng);
    test_avx2_kernels(rng);

    return 0;
}
