#include <iostream>
#include <immintrin.h>
#include <msgpack/fbuffer.hpp>
#include "assembled_chunk_msgpack.hpp"
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

uint64_t assembled_chunk::fpgacounts_begin() const {
    return isample * fpga_counts_per_sample;
}

uint64_t assembled_chunk::fpgacounts_end() const {
    return isample * this->fpga_counts_per_sample +
        constants::nt_per_assembled_chunk * this->fpga_counts_per_sample;
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

void assembled_chunk::downsample(assembled_chunk* dest,
                                 const assembled_chunk* src1,
                                 const assembled_chunk* src2) {

    float* dest_scales  = dest->scales;
    float* dest_offsets = dest->offsets;
    bool temp_scales = false;
    if (src1 == dest || src1 == dest) {
        // enable in-place downsampling
        dest_scales  = (float*)malloc(constants::nfreq_coarse_tot * nt_coarse * sizeof(float));
        dest_offsets = (float*)malloc(constants::nfreq_coarse_tot * nt_coarse * sizeof(float));
        temp_scales = true;
    }

    // Compute destination offset + scale.
    for (int i=0; i<(constants::nfreq_coarse_tot * nt_coarse); i++) {
        float scale_1  = src1->scales [i];
        float offset_1 = src1->offsets[i];
        float scale_2  = src2->scales [i];
        float offset_2 = src2->offsets[i];

        // Lower and upper limits of each offset, scale
        // choice... ignoring the fact that 0 and 255 are
        // marker values for masked values.
        float lo = min(offset_1, offset_2);
        float hi = max(offset_1 + scale_1 * 255., offset_2 + scale_2 * 255.);

        // Make the new range contain the old range.  (The
        // data may fill a smaller range than this...)
        dest_scales[i]  = (hi - lo) / 255.;
        dest_offsets[i] = lo;
    }

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float *scales_1  = src1->scales  + if_coarse * nt_coarse;
	const float *offsets_1 = src1->offsets + if_coarse * nt_coarse;
	const float *scales_2  = src2->scales  + if_coarse * nt_coarse;
	const float *offsets_2 = src2->offsets + if_coarse * nt_coarse;
	const float *scales_d  = dest_scales   + if_coarse * nt_coarse;
	const float *offsets_d = dest_offsets  + if_coarse * nt_coarse;

	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *data_1 = src1->data + if_fine * constants::nt_per_assembled_chunk;
	    const uint8_t *data_2 = src2->data + if_fine * constants::nt_per_assembled_chunk;
	    const uint8_t *data_d = dest->data + if_fine * constants::nt_per_assembled_chunk;

	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++) {
		float scale_1  = scales_1 [it_coarse];
		float offset_1 = offsets_1[it_coarse];
		float scale_2  = scales_2 [it_coarse];
		float offset_2 = offsets_2[it_coarse];

                float iscale_d = 1./scales_d[it_coarse];
                float offset_d = offsets_d[it_coarse];

		for (int it_fine = it_coarse*nt_per_packet; it_fine < (it_coarse+1)*nt_per_packet; it_fine++) {
		    float x1 = float(src_1[it_fine]);
		    float x2 = float(src_2[it_fine]);
                    float wt1 = ((x1==0) || (x1==255)) ? 0.0 : 1.0;
                    float wt2 = ((x2==0) || (x2==255)) ? 0.0 : 1.0;
                    float wtd = wt1 + wt2;
                    if (wtd == 0.)
                        src_d[it_fine] = 0;
                    else {
                        x1 = offset_1 + x1 * scale_1;
                        x2 = offset_2 + x2 * scale_2;
                        float xd = (x1*wt1 + x2*wt2) / wtd;
                        xd = (xd - offset_d) * iscale_d;
                        // FIXME -- round?
                        src_d[it_fine] = (uint8_t)lround(xd);
                    }
		}
	    }
	}
    }

    if (temp_scales) {
        memcpy(dest->scales,  dest_scales,  constants::nfreq_coarse_tot * nt_coarse * sizeof(float));
        memcpy(dest->offsets, dest_offsets, constants::nfreq_coarse_tot * nt_coarse * sizeof(float));
        free(dest_scales);
        free(dest_offsets);
    }
}


unique_ptr<assembled_chunk> assembled_chunk::make(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_, bool force_reference, bool force_fast)
{
    // FIXME -- if C++14 is available, use make_unique()
    if (force_reference && force_fast)
        throw runtime_error("ch_frb_io: assembled_chunk::make(): both force_reference and force_fast were set!");

#ifdef __AVX2__
    if (force_fast ||
        ((nt_per_packet_ == 16) && (nupfreq_ % 2 == 0) && !force_reference))
	return unique_ptr<fast_assembled_chunk>(new fast_assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_));
#else
    if (force_fast)
        throw runtime_error("ch_frb_io: assembled_chunk::make(): force_fast set on a machine without AVX2!");
#endif

    return unique_ptr<assembled_chunk>(new assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_));
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
    vector<hsize_t> scaleshape = { (hsize_t)constants::nfreq_coarse_tot,
                                   (hsize_t)this->nt_coarse };
    g_chunk.write_dataset("scales",  this->scales,  scaleshape);
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

void assembled_chunk::write_msgpack_file(const string &filename)
{
    FILE* f = fopen(filename.c_str(), "w+");
    if (!f)
        throw runtime_error("ch_frb_io: failed to open file " + filename + " for writing an assembled_chunk in msgpack format: " + strerror(errno));
    // msgpack buffer that will write to file "f"
    msgpack::fbuffer buffer(f);
    // Construct a shared_ptr from this, carefully
    shared_ptr<assembled_chunk> shthis(shared_ptr<assembled_chunk>(), this);
    msgpack::pack(buffer, shthis);
    if (fclose(f))
        throw runtime_error("ch_frb_io: failed to close assembled_chunk msgpack file " + filename + string(strerror(errno)));
}

shared_ptr<assembled_chunk> assembled_chunk::read_msgpack_file(const string &filename)
{
    struct stat st;
    if (stat(filename.c_str(), &st)) {
        throw runtime_error("ch_frb_io: failed to stat file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));
    }
    size_t len = st.st_size;
    FILE* f = fopen(filename.c_str(), "r");
    if (!f)
        throw runtime_error("ch_frb_io: failed to open file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));

    char* fdata = (char*)malloc(len);
    if (!fdata)
        throw runtime_error("ch_frb_io: failed to malloc an array of size " + to_string(len) + " for reading an assembled_chunk in msgpack format from file " + filename);

    size_t nr = fread(fdata, 1, len, f);
    if (nr != len)
        throw runtime_error("ch_frb_io: failed to read " + to_string(len) + " from file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));
    fclose(f);

    msgpack::object_handle oh = msgpack::unpack(fdata, len);
    msgpack::object obj = oh.get();
    shared_ptr<assembled_chunk> ch;
    //obj.convert(&ch);
    obj.convert(ch);
    free(fdata);
    return ch;
}


}  // namespace ch_frb_io
