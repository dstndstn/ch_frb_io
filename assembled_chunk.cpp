#include <iostream>
#include <immintrin.h>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


assembled_chunk::assembled_chunk(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_)
    : beam_id(beam_id_), 
      nupfreq(nupfreq_), 
      nt_per_packet(nt_per_packet_),
      fpga_counts_per_sample(fpga_counts_per_sample_), 
      nt_coarse(constants::nt_per_assembled_chunk / nt_per_packet),
      nscales(constants::nfreq_coarse_tot * nt_coarse),
      ndata(constants::nfreq_coarse_tot * nupfreq * constants::nt_per_assembled_chunk),
      ichunk(ichunk_),
      isample(ichunk * constants::nt_per_assembled_chunk)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("assembled_chunk constructor: bad beam_id argument");
    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("assembled_chunk constructor: bad nupfreq argument");
    if ((nt_per_packet <= 0) || !is_power_of_two(nt_per_packet) || (nt_per_packet > constants::nt_per_assembled_chunk))
	throw runtime_error("assembled_chunk constructor: bad nt_per_packet argument");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("assembled_chunk constructor: bad fpga_counts_per_sample argument");

    this->scales = aligned_alloc<float> (nscales);
    this->offsets = aligned_alloc<float> (nscales);
    this->data = aligned_alloc<uint8_t> (ndata);
}


assembled_chunk::~assembled_chunk()
{
    free(data);
    free(scales);
    free(offsets);
}


void assembled_chunk::fill_with_copy(const shared_ptr<assembled_chunk> &x)
{
    if (!x)
	throw runtime_error("assembled_chunk::fill_with_copy() called with empty pointer");
    if ((this->nupfreq != x->nupfreq) || (this->nt_per_packet != x->nt_per_packet))
	throw runtime_error("assembled_chunk::fill_with_copy() called on non-conformable chunks");

    if (x.get() == this)
	return;

    memcpy(this->data, x->data, ndata);
    memcpy(this->scales, x->scales, nscales * sizeof(float));
    memcpy(this->offsets, x->offsets, nscales * sizeof(float));
}


// Used in unit tests
void assembled_chunk::randomize(std::mt19937 &rng)
{
    for (int i = 0; i < ndata; i++) {
	// Assign ~10% probability to 0x00 or 0xff
	int x = randint(rng, -25, 281);
	x = max(x, 0);
	x = min(x, 255);
	this->data[i] = uint8_t(x);
    }

    uniform_rand(rng, this->scales, nscales);
    uniform_rand(rng, this->offsets, nscales);
}


// virtual member function; any changes made here should be reflected in override fast_assembled_chunk::add_packet().
void assembled_chunk::add_packet(const intensity_packet &packet)
{
    uint64_t packet_t0 = packet.fpga_count / uint64_t(fpga_counts_per_sample);

    // Offset relative to beginning of packet
    uint64_t t0 = packet_t0 - isample;
    
    // The runtime checks in intensity_network_stream::_process_packet() should
    // ensure that the following checks are redundant.  I decided to include the 
    // redundant checks here in the "generic" assembled_chunk::add_packet(), but 
    // omit them in fast_assembled_chunk::add_packet().

    bool bad = ((packet.nbeams != 1) ||
		(packet.nupfreq != this->nupfreq) ||
		(packet.ntsamp != this->nt_per_packet) ||
		(packet.fpga_counts_per_sample != this->fpga_counts_per_sample) ||
		(packet.fpga_count % (fpga_counts_per_sample * nt_per_packet)) ||
		(packet.beam_ids[0] != this->beam_id) ||
		(packet_t0 < isample) ||
		(packet_t0 + nt_per_packet > isample + constants::nt_per_assembled_chunk));

    if (_unlikely(bad))
	throw runtime_error("ch_frb_io: internal error in assembled_chunk::add_packet()");

    for (int f = 0; f < packet.nfreq_coarse; f++) {
	int coarse_freq_id = packet.coarse_freq_ids[f];

	this->scales[coarse_freq_id*nt_coarse + (t0/nt_per_packet)] = packet.scales[f];
	this->offsets[coarse_freq_id*nt_coarse + (t0/nt_per_packet)] = packet.offsets[f];

	for (int u = 0; u < nupfreq; u++) {
	    memcpy(data + (coarse_freq_id*nupfreq + u) * constants::nt_per_assembled_chunk + t0, 
		   packet.data + (f*nupfreq + u) * nt_per_packet,
		   nt_per_packet);
	}
    }
}


// virtual member function; any changes made here should be reflected in override fast_assembled_chunk::decode().
void assembled_chunk::decode(float *intensity, float *weights, int stride) const
{
    if (!intensity || !weights)
	throw runtime_error("ch_frb_io: null pointer passed to assembled_chunk::decode()");	
    if (stride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad stride passed to assembled_chunk::decode()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float *scales_f = this->scales + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;
	
	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * stride;
	    float *wt_f = weights + if_fine * stride;

	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++) {
		float scale = scales_f[it_coarse];
		float offset = offsets_f[it_coarse];
		
		for (int it_fine = it_coarse*nt_per_packet; it_fine < (it_coarse+1)*nt_per_packet; it_fine++) {
		    float x = float(src_f[it_fine]);
		    int_f[it_fine] = scale*x + offset;
		    wt_f[it_fine] = ((x==0) || (x==255)) ? 0.0 : 1.0;
		}
	    }
	}
    }
}


shared_ptr<assembled_chunk> assembled_chunk::make(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_)
{
#ifdef __AVX2__
    if ((nt_per_packet_ == 16) && (nupfreq_ % 2 == 0))
	return make_shared<fast_assembled_chunk> (beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_);
#endif

    return make_shared<assembled_chunk> (beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_);
}


void assembled_chunk::write_hdf5_file(const string &filename)
{
    bool write = true;
    bool clobber = true;
    hdf5_file f(filename, write, clobber);

    string chunkname = "/assembled-chunk-beam" + to_string(beam_id)
        + "-ichunk" + to_string(ichunk);
    bool create = true;
    hdf5_group g_chunk(f, chunkname, create);

    // Header
    g_chunk.write_attribute("beam_id", this->beam_id);
    g_chunk.write_attribute("nupfreq", this->nupfreq);
    g_chunk.write_attribute("nt_per_packet", this->nt_per_packet);
    g_chunk.write_attribute("fpga_counts_per_sample", this->fpga_counts_per_sample);
    g_chunk.write_attribute("nt_coarse", this->nt_coarse);
    g_chunk.write_attribute("nscales", this->nscales);
    g_chunk.write_attribute("ndata", this->ndata);
    g_chunk.write_attribute("ichunk", this->ichunk);
    g_chunk.write_attribute("isample", this->isample);

    // Offset & scale vectors
    vector<hsize_t> scaleshape = { (hsize_t)this->nscales };
    g_chunk.write_dataset("scales", this->scales, scaleshape);
    g_chunk.write_dataset("offsets", this->offsets, scaleshape);

    // Raw data
    int bitshuffle = 0;
    vector<hsize_t> datashape = {
        (hsize_t)constants::nfreq_coarse_tot,
        (hsize_t)nupfreq,
        (hsize_t)constants::nt_per_assembled_chunk };
    unique_ptr<hdf5_extendable_dataset<uint8_t> > data_dataset =
        make_unique<hdf5_extendable_dataset<uint8_t> >(g_chunk, "data", datashape, 2, bitshuffle);
    data_dataset->write(this->data, datashape);
    // close
    data_dataset = unique_ptr<hdf5_extendable_dataset<uint8_t> > ();
}



}  // namespace ch_frb_io
