#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


inline double intval(int ifreq, int ipol, int it)
{
    return 1.3183*ifreq + 1502.382*ipol + 0.328*it - 3201.91;
}

inline double wtval(int ifreq, int ipol, int it)
{
    return 1.021*ifreq - 1502.382*ipol - 0.283*it + 2038.18;
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    cout << "Warning: this will create a ~100 MB file in the current directory!\n"
	 << "If this is OK press return.  If not, press control-C!\n"
	 << "I AM WAITING, HUMAN: "
	 << flush;
    
    string dummy;
    getline(cin, dummy);

    const char *filename = "test_intensity_hdf5_file.hdf5";
    const int nfreq = 1024;
    const int npol = 2;
    const int nt_file = 16384;
    const double freq0 = 800.0;
    const double freq1 = 400.0;
    const double dt_sample = 1.0e-3;

    vector<string> pol = { "XX", "YY" };
    unique_ptr<intensity_hdf5_ofile> f = make_unique<intensity_hdf5_ofile> (filename, nfreq, pol, freq0, freq1, dt_sample);

    ssize_t expected_file_ipos = 0;
    vector<ssize_t> chunk_ipos_list;
    vector<ssize_t> chunk_nt_list;

    while (f->curr_nt < nt_file) {
	ssize_t chunk_ipos = expected_file_ipos + max(0, randint(-100,100));
	ssize_t chunk_nt = min((ssize_t)randint(1,33), nt_file - f->curr_nt);
	
	vector<float> intensity(nfreq * npol * chunk_nt);
	vector<float> weights(nfreq * npol * chunk_nt);

	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int ipol = 0; ipol < npol; ipol++) {
		for (int it = 0; it < chunk_nt; it++) {
		    int i = ifreq*npol*chunk_nt + ipol*chunk_nt + it;
		    intensity[i] = intval(ifreq, ipol, chunk_ipos + it);
		    weights[i] = wtval(ifreq, ipol, chunk_ipos + it);
		}
	    }
	}

	f->append_chunk(chunk_nt, &intensity[0], &weights[0], chunk_ipos);

	chunk_ipos_list.push_back(chunk_ipos);
	chunk_nt_list.push_back(chunk_nt);
	expected_file_ipos = chunk_ipos + chunk_nt;
    }

    return 0;
}
