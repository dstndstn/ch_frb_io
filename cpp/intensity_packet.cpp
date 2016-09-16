#include <cassert>
#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// Returns true if packet is good, false if bad.
// Note: no checking of freq_ids is performed!
bool intensity_packet::read(const uint8_t *src, int src_nbytes)
{
    if (_unlikely(src_nbytes < 24))
	return false;
    if (_unlikely(src_nbytes > constants::max_input_udp_packet_size))
	return false;

    memcpy(this, src, 24);

    if (_unlikely(protocol_version != 1))
	return false;
	
    // The following arithmetic operations have been chosen to avoid integer overflows,
    // and to minimize the number of branches required to validate the packet.
    uint64_t n1 = uint64_t(nbeams);
    uint64_t n2 = uint64_t(nfreq_coarse);
    uint64_t n3 = uint64_t(nupfreq);
    uint64_t n4 = uint64_t(ntsamp);

    // Expected header, data size
    uint64_t nh = 24 + 2*n1 + 2*n2 + 8*n1*n2;
    uint64_t nd = n1 * n2 * n3 * n4;

    if (_unlikely(uint64_t(src_nbytes) != nh+nd))
	return false;
    if (_unlikely(uint64_t(data_nbytes) != nd))
	return false;

    this->beam_ids = (uint16_t *) (src + 24);
    this->freq_ids = (uint16_t *) (src + 24 + 2*n1);
    this->scales = (float *) (src + 24 + 2*n1 + 2*n2);
    this->offsets = (float *) (src + 24 + 2*n1 + 2*n2 + 4*n1*n2);
    this->data = (uint8_t *) (src + nh);

    return true;
}


// No checking of the input is performed!
int intensity_packet::write(uint8_t *dst) const
{
    int n1 = this->nbeams;
    int n2 = this->nfreq_coarse;

    memcpy(dst, this, 24);
    memcpy(dst + 24, this->beam_ids, 2*n1);
    memcpy(dst + 24 + 2*n1, this->freq_ids, 2*n2);
    memcpy(dst + 24 + 2*n1 + 2*n2, this->scales, 4*n1*n2);
    memcpy(dst + 24 + 2*n1 + 2*n2 + 4*n1*n2, this->offsets, 4*n1*n2);
    memcpy(dst + 24 + 2*n1 + 2*n2 + 8*n1*n2, this->data, this->data_nbytes);

    return 24 + 2*n1 + 2*n2 + 8*n1*n2 + this->data_nbytes;
}


void test_packet_encoding()
{
    intensity_packet p;

    char *p0 = (char *) &p;
    assert((char *) &p.data_nbytes == p0+4);
    assert((char *) &p.fpga_counts_per_sample == p0+6);
    assert((char *) &p.fpga_count == p0+8);
    assert((char *) &p.nbeams == p0+16);
    assert((char *) &p.nfreq_coarse == p0+18);
    assert((char *) &p.nupfreq == p0+20);
    assert((char *) &p.ntsamp == p0+22);

    cerr << "test_packet_encoding(): success\n";
}


}   // namespace ch_frb_io
