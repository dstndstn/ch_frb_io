#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


assembled_chunk::assembled_chunk(int beam_id_, int nupfreq_, int fpga_counts_per_sample_, uint64_t chunk_t0_)
    : beam_id(beam_id), nupfreq(nupfreq_), fpga_counts_per_sample(fpga_counts_per_sample_), chunk_t0(chunk_t0_)
{
    if ((beam_id < 0) || (beam_id >= 65336))
	throw runtime_error("assembled_chunk constructor: bad beam_id argument");
    if ((nupfreq <= 0) || (nupfreq > 64))
	throw runtime_error("assembled_chunk constructor: bad nupfreq argument");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > 1024))
	throw runtime_error("assembled_chunk constructor: bad fpga_counts_per_sample argument");

    // FIXME hardcoded 1024 here
    this->intensity = aligned_alloc<float> (1024 * nfreq * nt_per_chunk);
    this->weights = aligned_alloc<float> (1024 * nfreq * nt_per_chunk);
}


~assembled_chunk::assembled_chunk()
{
    free(intensity);
    free(weights);
    this->intensity = this->weights = nullptr;
}



}  // namespace ch_frb_io
