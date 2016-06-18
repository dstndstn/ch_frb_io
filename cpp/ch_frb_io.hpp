#ifndef _CH_FRB_IO_HPP
#define _CH_FRB_IO_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++0x support (g++ -std=c++0x)"
#endif

#include <string>
#include <vector>

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif


struct noncopyable
{
    noncopyable() { }
    noncopyable(const noncopyable &) = delete;
    noncopyable& operator=(const noncopyable &) = delete;
};


struct intensity_hdf5_file : noncopyable {
    std::string filename;
    
    int nfreq;
    int npol;
    int nt_file;     // number of time samples in file (which can have gaps)
    int nt_logical;  // total number of samples in time range spanned by file, including gaps

    //
    // We currently throw an exception unless the frequencies are equally spaced and consecutive.
    // Therefore, we don't keep a full list of frequencies, just the endpoints of the range and number of samples.
    //
    // This still leaves two possibilities: the frequencies can either be ordered from lowest to highest,
    // or vice versa.  The 'frequencies_are_increasing' flag is set in the former case.  (Note that the
    // standard CHIME convention is frequencies_are_increasing=false.)
    //
    // The 'freq_lo_MHz' and 'freq_hi_MHz' fields are the lowest and highest frequencies in the entire
    // band.  Just to spell this out in detail, if frequencies_are_increasing=true, then the i-th channel
    // spans the frequency range
    //
    //   [ freq_lo + i*(freq_hi-freq_lo)/nfreq, freq_lo + (i+1)*(freq_hi-freq_lo)/nfreq ]
    //
    // and if frequences_are_increasing=false, then the i-th channel spans frequency range
    //
    //   [ freq_hi - (i+1)*(freq_hi-freq_lo)/nfreq, freq_hi - i*(freq_hi-freq_lo)/nfreq ]
    //
    // where 0 <= i <= nfreq-1.
    //
    bool frequencies_are_increasing;
    double freq_lo_MHz;
    double freq_hi_MHz;
    
    //
    // We distinguish between "file" time indices, which span the range 0 <= it < nt_file,
    // and "logical" time indices, which span the range 0 <= it < nt_logical.
    //
    // The i-th logical time sample spans the time range
    //
    //    [ time_lo + i*dt_sample, time_lo + (i+1)*dt_sample ]
    //
    double dt_sample;
    double time_lo;
    double time_hi;   // always equal to time_lo + nt_logical * dt_sample

    std::vector<double> times;                // 1D array of length nt_file
    std::vector<int> time_index_mapping;      // 1D array of length nt_file, which maps a "file" index to a "logical" index.

    // Polarization info (currently read from file but not really used)
    std::vector<std::string> polarizations;  // 1D array of length npol, each element is either "XX" or "YY"

    // 3d arrays of shape (nfreq, npol, nt_file), verbatim from the hdf5 file.
    // Rather than using them directly, you may want to use the member function get_unpolarized_intensity() below.
    std::vector<float> intensity;
    std::vector<float> weights;

    // Summary statistics
    double frac_ungapped;    // fraction of "logical" time samples which aren't in time gaps
    double frac_unmasked;    // fraction of _ungapped_ data with large weight

    // Construct from file.  If 'noisy' is true, then a one-line message will be printed when the file is read.
    explicit intensity_hdf5_file(const std::string &filename, bool noisy=true);

    //
    // Extracts a 2D array containing total intensity in time range [out_t0, out_t0+out_nt),
    // summing over polarizations.  The 'out_t0' and 'out_nt' are "logical" time indices, not
    // "file" time indices.  If this range of logical time indices contains gaps, the corresponding
    // entries of the 'out_int' and 'out_wt' arrays will be filled with zeros.
    // 
    // The 'out_int' and 'out_wt' arrays have shape (nfreq, out_nt).
    //
    // The 'out_stride' arg can be negative, if reversing the channel ordering is desired.
    // If out_stride is zero, it defaults to out_nt.
    //
    void get_unpolarized_intensity(float *out_int, float *out_wt, int out_t0, int out_nt, int out_stride=0) const;

    void run_unit_tests() const;
};


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_HPP
