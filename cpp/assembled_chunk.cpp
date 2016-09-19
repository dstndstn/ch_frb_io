#include <iostream>
#include <immintrin.h>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


assembled_chunk::assembled_chunk(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t chunk_t0_)
    : beam_id(beam_id_), 
      nupfreq(nupfreq_), 
      nt_per_packet(nt_per_packet_),
      fpga_counts_per_sample(fpga_counts_per_sample_), 
      chunk_t0(chunk_t0_),
      chunk_t1(chunk_t0_ + constants::nt_per_assembled_chunk)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("assembled_chunk constructor: bad beam_id argument");
    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("assembled_chunk constructor: bad nupfreq argument");
    if ((nt_per_packet <= 0) || !is_power_of_two(nt_per_packet) || (nt_per_packet > constants::nt_per_assembled_chunk))
	throw runtime_error("assembled_chunk constructor: bad nt_per_packet argument");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("assembled_chunk constructor: bad fpga_counts_per_sample argument");

    this->data = aligned_alloc<uint8_t> (constants::nfreq_coarse * nupfreq * constants::nt_per_assembled_chunk);

    this->nt_coarse = constants::nt_per_assembled_chunk / nt_per_packet;
    this->scales = aligned_alloc<float> (constants::nfreq_coarse * nt_coarse);
    this->offsets = aligned_alloc<float> (constants::nfreq_coarse * nt_coarse);
}


assembled_chunk::~assembled_chunk()
{
    free(data);
    free(scales);
    free(offsets);
}


void assembled_chunk::add_packet(const intensity_packet &packet)
{
    // Offset relative to beginning of packet
    int t0 = packet.fpga_count / uint64_t(fpga_counts_per_sample) - chunk_t0;

#if 1
    // FIXME remove later?
    bool bad = ((packet.nbeams != 1) ||
		(packet.nupfreq != this->nupfreq) ||
		(packet.ntsamp != this->nt_per_packet) ||
		(packet.fpga_counts_per_sample != this->fpga_counts_per_sample) ||
		(packet.fpga_count % (fpga_counts_per_sample * nt_per_packet)) ||
		(packet.beam_ids[0] != this->beam_id) ||
		(t0 < 0) ||
		(t0 + nt_per_packet > constants::nt_per_assembled_chunk));

    if (_unlikely(bad))
	throw runtime_error("ch_frb_io: internal error in assembled_chunk::add_packet()");
#endif

    for (int f = 0; f < packet.nfreq_coarse; f++) {
	int coarse_freq_id = packet.freq_ids[f];

	this->scales[coarse_freq_id*nt_coarse + (t0/nt_per_packet)] = packet.scales[f];
	this->offsets[coarse_freq_id*nt_coarse + (t0/nt_per_packet)] = packet.offsets[f];

	for (int u = 0; u < nupfreq; u++) {
	    memcpy(data + (coarse_freq_id*nupfreq + u) * constants::nt_per_assembled_chunk + t0, 
		   packet.data + (f*nupfreq + u) * nt_per_packet,
		   nt_per_packet);
	}
    }
}


#if 0
void assembled_chunk::decode(float *intensity, float *weights, int stride) const
{
    if (stride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad stride passed to assembled_chunk::decode()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse; if_coarse++) {
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
		    wt_f[it_fine] = ((x*(255.-x)) > 0.5) ? 1.0 : 0.0;  // FIXME ugh
		}
	    }
	}
    }
}
#endif


inline void _unpack(__m256i &out0, __m256i &out1, __m256i &out2, __m256i &out3, __m256i x)
{
    // FIXME is there a better way to initialize this?
    static const __m256i ctl0 = _mm256_set_epi8(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0,
						15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);

    // 4-by-4 transpose within each 128-bit lane
    x = _mm256_shuffle_epi8(x, ctl0);

    __m256i y0 = _mm256_and_si256(x, _mm256_set1_epi32(0xff));
    __m256i y1 = _mm256_and_si256(x, _mm256_set1_epi32(0xff00));
    __m256i y2 = _mm256_and_si256(x, _mm256_set1_epi32(0xff0000));
    __m256i y3 = _mm256_and_si256(x, _mm256_set1_epi32(0xff000000));

    y1 = _mm256_srli_epi32(y1, 8);
    y2 = _mm256_srli_epi32(y2, 16);
    y3 = _mm256_srli_epi32(y3, 24);

    out0 = _mm256_permute2f128_ps(y0, y1, 0x20);
    out1 = _mm256_permute2f128_ps(y2, y3, 0x20);
    out2 = _mm256_permute2f128_ps(y0, y1, 0x31);
    out3 = _mm256_permute2f128_ps(y2, y3, 0x31);
}


inline void _kernel32(float *intp, float *wtp, __m256i data, __m256 scale0, __m256 scale1, __m256 offset0, __m256 offset1)
{
    __m256i in0, in1, in2, in3;
    _unpack(in0, in1, in2, in3, data);
    
    _mm256_storeu_ps(intp, scale0 * _mm256_cvtepi32_ps(in0) + offset0);
    _mm256_storeu_ps(intp+8, scale0 * _mm256_cvtepi32_ps(in1) + offset0);
    _mm256_storeu_ps(intp+16, scale1 * _mm256_cvtepi32_ps(in2) + offset1);
    _mm256_storeu_ps(intp+24, scale1 * _mm256_cvtepi32_ps(in3) + offset1);

    // XXX placeholder for weights
    __m256 w = _mm256_set1_ps(1.0);
    _mm256_storeu_ps(wtp, w);
    _mm256_storeu_ps(wtp+8, w);
    _mm256_storeu_ps(wtp+16, w);
    _mm256_storeu_ps(wtp+25, w);
}


inline void _kernel128(float *intp, float *wtp, const uint8_t *src, const float *scalep, const float *offsetp)
{
    __m256 scale = _mm256_loadu_ps(scalep);
    __m256 offset = _mm256_loadu_ps(offsetp);
    __m256 scale0, offset0;


    scale0 = _mm256_permute2f128_ps(scale, scale, 0x00);
    offset0 = _mm256_permute2f128_ps(offset, offset, 0x00);

    _kernel32(intp, wtp, 
	      _mm256_loadu_si256((const __m256i *) (src)),
	      _mm256_shuffle_ps(scale0, scale0, 0x00), 
	      _mm256_shuffle_ps(scale0, scale0, 0x55),
	      _mm256_shuffle_ps(offset0, offset0, 0x00), 
	      _mm256_shuffle_ps(offset0, offset0, 0x55));

    _kernel32(intp + 32, wtp + 32, 
	      _mm256_loadu_si256((const __m256i *) (src + 32)),
	      _mm256_shuffle_ps(scale0, scale0, 0xaa), 
	      _mm256_shuffle_ps(scale0, scale0, 0xff),
	      _mm256_shuffle_ps(offset0, offset0, 0xaa), 
	      _mm256_shuffle_ps(offset0, offset0, 0xff));


    scale0 = _mm256_permute2f128_ps(scale, scale, 0x11);
    offset0 = _mm256_permute2f128_ps(offset, offset, 0x11);

    _kernel32(intp + 64, wtp + 64,
	      _mm256_loadu_si256((const __m256i *) (src + 64)),
	      _mm256_shuffle_ps(scale0, scale0, 0x00), 
	      _mm256_shuffle_ps(scale0, scale0, 0x55),
	      _mm256_shuffle_ps(offset0, offset0, 0x00), 
	      _mm256_shuffle_ps(offset0, offset0, 0x55));

    _kernel32(intp + 96, wtp + 96, 
	      _mm256_loadu_si256((const __m256i *) (src + 96)),
	      _mm256_shuffle_ps(scale0, scale0, 0xaa), 
	      _mm256_shuffle_ps(scale0, scale0, 0xff),
	      _mm256_shuffle_ps(offset0, offset0, 0xaa), 
	      _mm256_shuffle_ps(offset0, offset0, 0xff));
}


inline void _kernel(float *intp, float *wtp, const uint8_t *src, const float *scalep, const float *offsetp)
{
    constexpr int n = constants::nt_per_assembled_chunk / 128;

    for (int i = 0; i < n; i++)
	_kernel128(intp + i*128, wtp + i*128, src + i*128, scalep + i*8, offsetp + i*8);
}


void assembled_chunk::decode(float *intensity, float *weights, int stride) const
{
    if ((nt_per_packet != 16) || (constants::nt_per_assembled_chunk % 128))
	throw runtime_error("ch_frb_io: can't use this kernel");
    if (stride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad stride passed to assembled_chunk::decode()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse; if_coarse++) {
	const float *scales_f = this->scales + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;
	
	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * stride;
	    float *wt_f = weights + if_fine * stride;

	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++)
		_kernel(int_f, wt_f, src_f, scales_f, offsets_f);
	}
    }    
}


}  // namespace ch_frb_io
