#include <iomanip>
#include <immintrin.h>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// Utils


template<unsigned int N>
inline float _extract128(__m128 x)
{
    // _mm_extract_ps() returns int instead of float?!
    union { int i; float x; } u;
    u.i = _mm_extract_ps(x, N);
    return u.x;
}

template<unsigned int N>
inline float _extract256(__m256 x)
{
    __m128 x2 = _mm256_extractf128_ps(x, N/4);
    return _extract128<N%4> (x2);
}

template<unsigned int N> inline void _vstr8_partial(stringstream &ss, __m256i x, bool hexflag);
template<unsigned int N> inline void _vstr32_partial(stringstream &ss, __m256i x, bool hexflag);
template<unsigned int N> inline void _vstr_partial(stringstream &ss, __m256 x);

template<> inline void _vstr8_partial<0>(stringstream &ss, __m256i x, bool hex) { return; }
template<> inline void _vstr32_partial<0>(stringstream &ss, __m256i x, bool hex) { return; }
template<> inline void _vstr_partial<0>(stringstream &ss, __m256 x) { return; }

template<unsigned int N> 
inline void _vstr8_partial(stringstream &ss, __m256i x, bool hexflag) 
{
    _vstr8_partial<N-1>(ss, x, hexflag);
    if (hexflag)
	ss << " " << setfill('0') << setw(2) << hex << uint32_t(uint8_t(_mm256_extract_epi8(x,N-1)));
    else
	ss << " " << int32_t(_mm256_extract_epi8(x,N-1));
}

template<unsigned int N>
inline void _vstr32_partial(stringstream &ss, __m256i x, bool hexflag) 
{
    _vstr32_partial<N-1>(ss, x, hexflag);
    if (hexflag)
	ss << " " << setfill('0') << setw(8) << hex << uint32_t(_mm256_extract_epi32(x,N-1));
    else
	ss << " " << _mm256_extract_epi32(x,N-1);
}

template<unsigned int N>
inline void _vstr_partial(stringstream &ss, __m256 x)
{
    _vstr_partial<N-1>(ss, x);
    ss << " " << _extract256<N-1>(x);
}


inline string _vstr8(__m256i x, bool hexflag=true)
{
    stringstream ss;
    ss << "[";
    _vstr8_partial<32> (ss, x, hexflag);
    ss << " ]";
    return ss.str();
}

inline string _vstr32(__m256i x, bool hexflag=false)
{
    stringstream ss;
    ss << "[";
    _vstr32_partial<8> (ss, x, hexflag);
    ss << " ]";
    return ss.str();
}

inline string _vstr(__m256 x)
{
    stringstream ss;
    ss << "[";
    _vstr_partial<8> (ss, x);
    ss << " ]";
    return ss.str();
}


// -------------------------------------------------------------------------------------------------
//
// add_packet_kernel


inline void _add_packet_kernel(uint8_t *dst, const uint8_t *src, int nupfreq)
{
    constexpr int s = constants::nt_per_assembled_chunk;

    for (int i = 0; i < nupfreq; i += 2) {
	__m256i x = _mm256_loadu_si256((const __m256i *) (src + 16*i));
	__m128i x0 = _mm256_extractf128_si256(x, 0);
	__m128i x1 = _mm256_extractf128_si256(x, 1);
	
	_mm_storeu_si128((__m128i *) (dst + i*s), x0);
	_mm_storeu_si128((__m128i *) (dst + (i+1)*s), x1);
    }
}


// -------------------------------------------------------------------------------------------------
//
// Decode kernels


inline void _decode_unpack(__m256i &out0, __m256i &out1, __m256i &out2, __m256i &out3, __m256i x)
{
    // FIXME is there a better way to initialize this?
    static const __m256i ctl0 = _mm256_set_epi8(15,11,7,3,14,10,2,6,13,9,5,1,12,8,4,0,
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

    out0 = _mm256_permute2f128_si256(y0, y1, 0x20);
    out1 = _mm256_permute2f128_si256(y2, y3, 0x20);
    out2 = _mm256_permute2f128_si256(y0, y1, 0x31);
    out3 = _mm256_permute2f128_si256(y2, y3, 0x31);
}


inline void _decode_weights(float *wtp, __m256i x, __m256i i0, __m256i i254, __m256 f0, __m256 f1)
{
    __m256i gt0 = _mm256_cmpgt_epi32(x, i0);
    __m256i gt254 = _mm256_cmpgt_epi32(x, i254);
    __m256i valid = _mm256_andnot_si256(gt254, gt0);
    __m256  wt = _mm256_blendv_ps(f0, f1, (__m256)valid);

    _mm256_storeu_ps(wtp, wt);
}


inline void _decode_kernel32(float *intp, float *wtp, __m256i data, __m256 scale0, __m256 scale1, __m256 offset0, __m256 offset1)
{
    __m256i in0, in1, in2, in3;
    _decode_unpack(in0, in1, in2, in3, data);
    
    _mm256_storeu_ps(intp, scale0 * _mm256_cvtepi32_ps(in0) + offset0);
    _mm256_storeu_ps(intp+8, scale0 * _mm256_cvtepi32_ps(in1) + offset0);
    _mm256_storeu_ps(intp+16, scale1 * _mm256_cvtepi32_ps(in2) + offset1);
    _mm256_storeu_ps(intp+24, scale1 * _mm256_cvtepi32_ps(in3) + offset1);

    __m256i i0 = _mm256_set1_epi32(0);
    __m256i i254 = _mm256_set1_epi32(254);
    __m256 f0 = _mm256_set1_ps(0.0);
    __m256 f1 = _mm256_set1_ps(1.0);

    _decode_weights(wtp, in0, i0, i254, f0, f1);
    _decode_weights(wtp+8, in1, i0, i254, f0, f1);
    _decode_weights(wtp+16, in2, i0, i254, f0, f1);
    _decode_weights(wtp+24, in3, i0, i254, f0, f1);
}


inline void _decode_kernel128(float *intp, float *wtp, const uint8_t *src, const float *scalep, const float *offsetp)
{
    __m256 scale = _mm256_loadu_ps(scalep);
    __m256 offset = _mm256_loadu_ps(offsetp);
    __m256 scale0, offset0;

    scale0 = _mm256_permute2f128_ps(scale, scale, 0x00);
    offset0 = _mm256_permute2f128_ps(offset, offset, 0x00);

    _decode_kernel32(intp, wtp, 
		     _mm256_loadu_si256((const __m256i *) (src)),
		     _mm256_shuffle_ps(scale0, scale0, 0x00), 
		     _mm256_shuffle_ps(scale0, scale0, 0x55),
		     _mm256_shuffle_ps(offset0, offset0, 0x00), 
		     _mm256_shuffle_ps(offset0, offset0, 0x55));

    _decode_kernel32(intp + 32, wtp + 32, 
		     _mm256_loadu_si256((const __m256i *) (src + 32)),
		     _mm256_shuffle_ps(scale0, scale0, 0xaa), 
		     _mm256_shuffle_ps(scale0, scale0, 0xff),
		     _mm256_shuffle_ps(offset0, offset0, 0xaa), 
		     _mm256_shuffle_ps(offset0, offset0, 0xff));


    scale0 = _mm256_permute2f128_ps(scale, scale, 0x11);
    offset0 = _mm256_permute2f128_ps(offset, offset, 0x11);

    _decode_kernel32(intp + 64, wtp + 64,
		     _mm256_loadu_si256((const __m256i *) (src + 64)),
		     _mm256_shuffle_ps(scale0, scale0, 0x00), 
		     _mm256_shuffle_ps(scale0, scale0, 0x55),
		     _mm256_shuffle_ps(offset0, offset0, 0x00), 
		     _mm256_shuffle_ps(offset0, offset0, 0x55));

    _decode_kernel32(intp + 96, wtp + 96, 
		     _mm256_loadu_si256((const __m256i *) (src + 96)),
		     _mm256_shuffle_ps(scale0, scale0, 0xaa), 
		     _mm256_shuffle_ps(scale0, scale0, 0xff),
		     _mm256_shuffle_ps(offset0, offset0, 0xaa), 
		     _mm256_shuffle_ps(offset0, offset0, 0xff));
}


inline void _decode_kernel(float *intp, float *wtp, const uint8_t *src, const float *scalep, const float *offsetp)
{
    static_assert(constants::nt_per_assembled_chunk % 128 == 0, "_decode_kernel() assumes nt_per_assembled_chunk divisible by 128");

    constexpr int n = constants::nt_per_assembled_chunk / 128;

    for (int i = 0; i < n; i++)
	_decode_kernel128(intp + i*128, wtp + i*128, src + i*128, scalep + i*8, offsetp + i*8);
}


// -------------------------------------------------------------------------------------------------
//
// class fast_assembled_chunk


fast_assembled_chunk::fast_assembled_chunk(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t chunk_t0_) :
    assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, chunk_t0_)
{
    if (nt_per_packet_ != 16)
	throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk constructor called with nt_per_packet != 16");
    if (nupfreq % 2 != 0)
	throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk constructor called with odd value of nupfreq");
}


// virtual override
void fast_assembled_chunk::add_packet(const intensity_packet &packet)
{
    // Offset relative to beginning of packet
    int t0 = packet.fpga_count / uint64_t(fpga_counts_per_sample) - chunk_t0;

    for (int f = 0; f < packet.nfreq_coarse; f++) {
	int coarse_freq_id = packet.freq_ids[f];

	int d = coarse_freq_id*nt_coarse + (t0/nt_per_packet);
	this->scales[d] = packet.scales[f];
	this->offsets[d] = packet.offsets[f];

	uint8_t *dst = data + coarse_freq_id * nupfreq * constants::nt_per_assembled_chunk + t0;
	const uint8_t *src = packet.data + f * nupfreq * 16;

	_add_packet_kernel(dst, src, nupfreq);
    }
}


// virtual override
void fast_assembled_chunk::decode(float *intensity, float *weights, int stride) const
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

	    _decode_kernel(int_f, wt_f, src_f, scales_f, offsets_f);
	}
    }    
}


// -------------------------------------------------------------------------------------------------
//
// Testing


// helper function used by test_fast_decode_kernel()
static vector<float> randvec(std::mt19937 &rng, ssize_t n)
{
    vector<float> ret(n);
    uniform_rand(rng, &ret[0], n);
    return ret;
}


void peek_at_unpack_kernel()
{
    uint8_t v[32];
    for (int i = 0; i < 32; i++)
	v[i] = i+1;

    __m256i x = _mm256_loadu_si256((const __m256i *) v);

    __m256i y0, y1, y2, y3;
    _decode_unpack(y0, y1, y2, y3, x);

    cout << _vstr8(x,false) << endl
	 << _vstr32(y0,false) << endl
	 << _vstr32(y1,false) << endl
	 << _vstr32(y2,false) << endl
	 << _vstr32(y3,false) << endl;
}


void test_fast_decode_kernel(std::mt19937 &rng)
{
    cerr << "test_fast_decode_kernel()";

    // Required by fast decode kernel
    const int nt_per_packet = 16;

    // Initialized arbitrarily, since decode() doesn't use them.
    const int beam_id = 0;
    const int fpga_counts_per_sample = 384;
    const uint64_t chunk_t0 = 0;

    for (int iouter = 0; iouter < 128; iouter++) {
	cerr << ".";

	// Randomized in every iteration
	const int nupfreq = randint(rng, 1, 17);
	const int stride = randint(rng, constants::nt_per_assembled_chunk, constants::nt_per_assembled_chunk + 16);
	
	auto chunk0 = make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);
	auto chunk1 = make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, chunk_t0);

	chunk0->randomize(rng);
	chunk1->fill_with_copy(chunk0);

	int nfreq_fine = constants::nfreq_coarse * nupfreq;
	vector<float> intensity0 = randvec(rng, nfreq_fine * stride);
	vector<float> intensity1 = randvec(rng, nfreq_fine * stride);
	vector<float> weights0 = randvec(rng, nfreq_fine * stride);
	vector<float> weights1 = randvec(rng, nfreq_fine * stride);

	chunk0->decode(&intensity0[0], &weights0[0], stride);
	chunk1->decode(&intensity1[0], &weights1[0], stride);

	for (int ifreq = 0; ifreq < nfreq_fine; ifreq++) {
	    for (int it = 0; it < constants::nt_per_assembled_chunk; it++) {
		int i = ifreq*stride + it;
		int j = ifreq*constants::nt_per_assembled_chunk + it;

		if (fabs(intensity0[i] - intensity1[i]) > 1.0e-5) {
		    cerr << "\n " << nupfreq << " " << ifreq << " " << it 
			 << " " << int32_t(chunk0->data[j]) << " " << int32_t(chunk1->data[j])
			 << " " << intensity0[i] << " " << intensity1[i] << endl;
		    throw runtime_error("test_fast_decode_kernel: intensity mismatch");
		}

		if (weights0[i] != weights1[i]) {
		    cerr << "\n " << nupfreq << " " << ifreq << " " << it 
			 << " " << int32_t(chunk0->data[j]) << " " << int32_t(chunk1->data[j])
			 << " " << weights0[i] << " " << weights1[i] << endl;
		    throw runtime_error("test_fast_decode_kernel: weights mismatch");
		}
	    }
	}
    }

    cerr << "success\n";
}


}  // namespace ch_frb_io
