//
// The following pseudo header file isn't actually valid C++ code, but is
// intended to document the L0_L1 UDP packet format in a C++-like syntax.
// For real packet-encoding code, see intensity_network_ostream.cpp.
//
// We don't bother network-encoding ints and floats, under the assumption
// that all nodes in the network are Intel CPUs with little-endian ints
// and IEEE 754 floats.
//
// The packet header size in bytes is currently
//    24 + 2*nbeam + 2*nfreq + 8*nbeam*nfreq
//
// For full CHIME, we expect to use (nbeam, nfreq, nupfreq, ntsamp) = (8, 4, 16, 16)
// which gives a 304-byte header and an 8192-byte data segment.
//
// For the pathfinder, I'm currently using  (nbeam, nfreq, nupfreq, ntsamp) = (1, 32, 1, 256) 
// which gives a 346-byte header and an 8192-byte data segment.
//

struct L0_L1_header {
    // This file describes protocol version 1.
    // A 32-bit protocol number is overkill, but means that all fields below are aligned on their
    // "natural" boundaries (i.e. fields with size Nbytes have byte offsets which are multiples of Nbytes)
    uint32_t    protocol_version;    // byte offset 0
    
    // Size of the 'data' array in bytes, where a negative size means "bitshuffle-compressed".
    // If positive, must equal (nbeam * nfreq * nupfreq * ntsamp)
    int16_t     data_nbytes;        // byte offset 4

    // This is the duration of a time sample, in FPGA counts.
    // The duration in seconds is dt = (2.56e-6 * fpga_counts_per_sample)
    uint16_t    fpga_counts_per_sample;   // byte offset 6

    // This is the time index (in FPGA counts) of the first time sample in the packet.
    // The packet sender is responsible for "unwrapping" the 32-bit FPGA timestamp to 64 bits.
    // Must be divisible by fpga_counts_per_sample.
    uint64_t    fpga_count;        // byte offset 8

    uint16_t    nbeam;      // number of beams in packet, byte offset 16
    uint16_t    nfreq;      // number of "FPGA" frequency channels between 0 and 1024, byte offset 18
    uint16_t    nupfreq;    // upsampling factor (probably either 1 or 16), byte offset 20
    uint16_t    ntsamp;     // number of time samples in packet, byte offset 22

    // The beams and FPGA freqencies in a particular (L0,L1) pair are not assumed contiguous.
    uint16_t    beam_ids[nbeam];   // beam_ids are between 0 and 1024 for full chime, byte offset of array is 24
    uint16_t    freq_ids[nfreq];   // freq_ids are between 0 and 1024 for full chime, byte offset of array is (24+2*nbeam)
    
    //
    // The intensities are 8-bit encoded with 'scale' and 'offset' parameters:
    //    Floating-point intensity = scale * (8-bit value) + offset
    //
    // A design decision: what granularity should the scale and offset parameters have?
    // I thought it made sense to allow the scale and offset to depend on both beam and
    // frequency.  This increases the header size so that the header is ~4% of the total
    // packet.  While not negligible, this increase seemed small enough to be a reasonable
    // price for the robustness of per-(beam,frequency) encoding.
    //
    float32     scale[nbeam * nfreq];   // byte offset (24 + 2*nbeam + 2*nfreq)
    float32     offset[nbeam * nfreq];  // byte offset (24 + 2*nbeam + 2*nfreq + 2*nbeam*nfreq)

    // 8-bit encoded data, as a 4D array with shape (nbeam, nfreq, nupfreq, ntsamp).
    // We need a sentinel value to indicate "this entry in the array is invalid".
    // We currently take both endpoint values to be sentinels, i.e. if a byte is either
    // 0x00 or 0xff, it should be treated as masked.
    uint8_t     data[ abs(data_nbytes) ];   // byte offset (24 + 2*nbeam + 2*nfreq + 4*nbeam*nfreq)
};
