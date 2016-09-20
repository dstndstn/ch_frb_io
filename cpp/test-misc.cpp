#include "ch_frb_io_internals.hpp"

using namespace ch_frb_io;

int main(int argc, char **argv)
{
    test_lexical_cast();
    test_packet_encoding();

    // peek_at_avx2_kernels();
    return 0;
}
