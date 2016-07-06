#include <memory>
#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


intensity_hdf5_ofile::intensity_hdf5_ofile(const string &filename_, int nfreq_, const vector<string> &pol, double freq0, double freq1,
					   double dt_sample_, ssize_t ipos0, double time0, int bitshuffle, int nt_chunk)
    : filename(filename_), 
      dt_sample(dt_sample_),
      nfreq(nfreq_), 
      npol(pol.size()),
      curr_nt(0),
      curr_ipos(ipos0),
      curr_time(time0)
{
    // A lot of argument checking

    if (nfreq <= 0)
	throw runtime_error(filename + ": expected nfreq > 0 in intensity_hdf5_ofile constructor");
    if (npol == 0)
	throw runtime_error(filename + ": expected pol.size() > 0 in intensity_hdf5_ofile constructor");
    if (dt_sample <= 0.0)
	throw runtime_error(filename + ": expected dt_sample > 0 in intensity_hdf5_ofile constructor");
    if ((bitshuffle < 0) || (bitshuffle > 3))
	throw runtime_error(filename + ": expected bitshuffle = 0,1,2 or 3 in intensity_hdf5_ofile constructor");
    if (nt_chunk <= 0)
	throw runtime_error(filename + ": expected nt_chunk > 0 in intensity_hdf5_ofile constructor");

    for (int ipol = 0; ipol < npol; ipol++) {
	if ((pol[ipol] != "XX") && (pol[ipol] != "YY"))
	    throw runtime_error(filename + ": expected all polarization strings to be either 'XX' or 'YY' in intensity_hdf5_ofile constructor");
	
	for (int jpol = 0; jpol < ipol; jpol++)
	    if (pol[ipol] == pol[jpol])
		throw runtime_error(filename + ": duplicate polarizations in intensity_hdf5_ofile constructor");
    }

    bool write = true;
    bool clobber = true;
    hdf5_file f(filename, write, clobber);

    bool create = true;
    hdf5_group g_index_map(f, "index_map", create);
    hdf5_group g_root(f, ".");

    vector<double> freq(nfreq);
    for (int ifreq = 0; ifreq < nfreq; ifreq++)
	freq[ifreq] = (ifreq*freq1 + (nfreq-ifreq)*freq0) / (double)nfreq;

    vector<hsize_t> freq_shape = { (hsize_t)nfreq };
    g_index_map.write_dataset("freq", &freq[0], freq_shape);
    
    vector<hsize_t> pol_shape = { pol.size() };
    g_index_map.write_string_dataset("pol", pol, pol_shape);
    
    // No bitshuffle compression
    vector<hsize_t> tchunk_shape = { 16384 };
    this->time_dataset = make_unique<hdf5_extendable_dataset<double> >(g_index_map, "time", tchunk_shape, 0);

    vector<hsize_t> chunk_shape = { (hsize_t)nfreq, 2, (hsize_t)nt_chunk };
    this->intensity_dataset = make_unique<hdf5_extendable_dataset<float> > (g_root, "intensity", chunk_shape, 2, bitshuffle);
    this->weights_dataset = make_unique<hdf5_extendable_dataset<float> > (g_root, "weight", chunk_shape, 2, bitshuffle);
}


void intensity_hdf5_ofile::append_chunk(ssize_t nt_chunk, float *intensity, float *weights, ssize_t chunk_ipos, double chunk_t0)
{
    // Argument checking
    
    if (nt_chunk <= 0)
	throw runtime_error(filename + ": append_chunk() was called with nt_chunk <= 0");
    if (!intensity || !weights)
	throw runtime_error(filename + ": append_chunk() was called with null pointer");

    if (this->curr_nt > 0) {
	if (chunk_ipos < this->curr_ipos)
	    throw runtime_error(filename + ": append_chunk() was called with non-increasing or overlapping sample_ipos ranges");

	double expected_chunk_t0 = curr_time + dt_sample * (chunk_ipos - curr_ipos);
	double allowed_drift = 1.0e-2 * dt_sample;
	
	if (fabs(chunk_t0 - expected_chunk_t0) > allowed_drift)
	    throw runtime_error(filename + ": append_chunk() was called with wrong timestamp or too much timestamp drift");
    }

    // Write data

    vector<double> times(nt_chunk);
    for (int it = 0; it < nt_chunk; it++)
	times[it] = chunk_t0 + it * dt_sample;

    vector<hsize_t> tchunk_shape = { (hsize_t) nt_chunk };
    this->time_dataset->write(&times[0], tchunk_shape);

    vector<hsize_t> chunk_shape = { (hsize_t)nfreq, (hsize_t)npol, (hsize_t)nt_chunk };
    this->intensity_dataset->write(intensity, chunk_shape);
    this->weights_dataset->write(weights, chunk_shape);

    // Update state

    this->curr_nt += nt_chunk;
    this->curr_ipos = chunk_ipos + nt_chunk;
    this->curr_time = chunk_t0 + nt_chunk * dt_sample;
}


void intensity_hdf5_ofile::append_chunk(ssize_t nt_chunk, float *intensity, float *weight, ssize_t chunk_ipos)
{
    double chunk_t0 = this->curr_time + dt_sample * (chunk_ipos - curr_ipos);
    this->append_chunk(nt_chunk, intensity, weight, chunk_ipos, chunk_t0);
}


}  // namespace ch_frb_io
