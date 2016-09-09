#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


static void usage()
{
    cerr << "usage: ch-show-intensity-file [-t] [-m] <filename1.h5> <filename2.h5> ...\n"
	 << "       if -t is specified, then unit tests will be run\n"
	 << "       if -m is specified, then the mean intensity/weight will be shown in each channel (generates a lot of output!)\n";

    exit(2);
}


int main(int argc, char **argv)
{
    vector<string> filename_list;
    bool run_unit_tests = false;
    bool show_mean_values = false;

    // Low-budget command line parsing

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t"))
	    run_unit_tests = true;
        else if (!strcmp(argv[i], "-m"))
	    show_mean_values = true;
	else
	    filename_list.push_back(string(argv[i]));
    }

    if (filename_list.size() == 0)
	usage();

    for (unsigned int ifile = 0; ifile < filename_list.size(); ifile++) {
	if (ifile > 1)
	    cout << "\n";

	intensity_hdf5_file f(filename_list[ifile], true);  // noisy=true
    
	if (show_mean_values) {
	    int ns = f.npol * f.nt_file;

	    if (ns == 0)
		throw runtime_error("No data in file!\n");

	    for (int ifreq = 0; ifreq < f.nfreq; ifreq++) {
		double sum_wi = 0.0;
		double sum_wt = 0.0;

		for (int s = 0; s < ns; s++) {
		    sum_wi += f.weights[ifreq*ns+s] * f.intensity[ifreq*ns+s];
		    sum_wt += f.weights[ifreq*ns+s];
		}

		cout << "channel " << ifreq << "/" << f.nfreq << ": ";
		
		if (sum_wt <= 0.0) {
		    cout << "empty\n";
		    continue;
		}

		cout << " mean_int=" << (sum_wi/sum_wt) << ", mean_wt=" << (sum_wt/ns) << endl;
	    }
	}

	cout << "    nfreq = " << f.nfreq << endl
	     << "    npol = " << f.npol << endl
	     << "    nt_file = " << f.nt_file << "                  # Number of time indices in file\n"
	     << "    nt_logical = " << f.nt_logical << "               # Number of \"logical\" time indices spanned by file, including missing time regions\n"
	     << "    frequencies_are_increasing = " << f.frequencies_are_increasing << "   # boolean (should be 'false' for CHIME)\n" 
	     << "    freq_lo_MHz = " << f.freq_lo_MHz << endl
	     << "    freq_hi_MHz = " << f.freq_hi_MHz << endl
	     << "    time_lo = " << f.time_lo << endl
	     << "    time_hi = " << f.time_hi << endl
	     << "    dt_sample = " << f.dt_sample << endl
	     << "    frac_ungapped = " << f.frac_ungapped << "     # Fraction of time interval which is present in file (i.e. nt_file/nt_logical)\n"
	     << "    frac_unmasked = " << f.frac_unmasked << "     # This is sum(weights)/max(weights) with time gaps omitted\n";
	
	if (run_unit_tests)
	    f.run_unit_tests();
    }

    return 0;
}
