#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


static void usage()
{
    cerr << "usage: ch-show-intensity-file [-t] <filename1.h5> <filename2.h5> ...\n"
	 << "       if the -t flag is specified, then unit tests will be run\n";

    exit(2);
}


int main(int argc, char **argv)
{
    vector<string> filename_list;
    bool run_unit_tests = false;

    // Low-budget command line parsing

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t"))
	    run_unit_tests = true;
	else
	    filename_list.push_back(string(argv[i]));
    }

    if (filename_list.size() == 0)
	usage();

    for (unsigned int i = 0; i < filename_list.size(); i++) {
	if (i > 1)
	    cout << "\n";

	bool noisy = true;
	intensity_hdf5_file f(filename_list[i], noisy);
    
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
