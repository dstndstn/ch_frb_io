ch_frb_io: C++/python library for CHIME-FRB file and network streams.

CHIME FRB data has a similar format to CHIME cosmology data and can be read
using the tools in 'caput', in particular the caput.tod module
(https://github.com/radiocosmology/caput).
This package contains additional tools, such as simple file data streamers and a C
interface.  

There are now some minimal C++ classes for network input/output streams,
with a few loose ends which should be improved.

Currently there is no shared code between the C++ and python parts, so this is really two
independent libraries in the same git repository.  In fact the C++ and python parts use
different build systems, so you have to build and install them independently.  This is
something that should probably be fixed later!

This code is mostly used as a library in other places (e.g. https://github.com/kmsmith137/ch_vdif_assembler
or https://github.com/kmsmith137/rf_pipelines), but there are a few standalone command-line utilities here:
```
ch-show-intensity-file   print summary statistics for an hdf5 intensity file
ch-plot-intensity-file   make a quick waterfall plot from an hdf5 intensity file (requires rf_pipelines)
decompress-chfrb-data    bitshuffle-decompress an hdf5 intensity file
```


### DEPENDENCIES

  1. libhdf5 (https://www.hdfgroup.org/HDF5/release/obtain5.html)
     Note that this is a link to HDF5 v1.8.  I imagine v1.10 also works but haven't tested it yet.

  2. Optional but recommended: bitshuffle (https://github.com/kiyo-masui/bitshuffle)
     You'll need this if you want to use bitshuffle-compressed files (note that CHIME pathfinder
     data is generally bitshuffle-compresed).

     A hint for installing bitshuffle:
     ```
     git clone https://github.com/kiyo-masui/bitshuffle.git
     cd bitshuffle/

     # The HDF5 library can dynamically load the bitshuffle plugin, i.e. you don't need
     # to link the bitshuffle library when you compile ch_frb_io, but you need to set this
     # environment variable to tell libhdf5 where to look.  Suggest adding this to .bashrc!

     export HDF5_PLUGIN_PATH=$HOME/lib/hdf5_plugins

     # If you have root privs and want to install "system-wide", omit the --user flag
     # The --h5plugin* flags will build/install the plugin needed to use bitshuffle from C++

     python setup.py install --user --h5plugin --h5plugin-dir=$HOME/lib/hdf5_plugins
     ```


### INSTALLATION (C++)

  - `cd cpp/`

  - Create a file ./Makefile.local containing compiler flags, library locations, etc.
    The details are described in the Makefile.  There are some examples in the site/
    directory.  (You may be able to just symlink one of these examples to ./Makefile.local)

  - `make all install`

  - If you want to run the unit test (warning: creates 100MB temp file in current directory):
    `./test-intensity-hdf5-file`


INSTALLATION (PYTHON)
---------------------

  - To build and install, do: `python setup.py install --user`
    (If you have root privs and want to install "system-wide", omit the --user flag)


### LOOSE ENDS IN NETWORKING CODE

  - Currently, we have to run `test-network-streams.cpp` at very low throughput (0.1 Gbps)
    to avoid dropping packets.  This means that the unit tests take about an hour to run,
    which isn't really a problem, but is indicative of deeper performance problems?  It
    would be nice to understand where the bottleneck is.

  - We don't have thread-safe logging, so diagnostic messages from the various threads
    sometimes interleave each other and are unreadable.  (Implementing thread-safe
    logging is a general todo item for the whole project, not just libch_frb_io!)

  - Cleanup: it would be great to switch from pthreads to C++11 threads.  (The C++11 API
    is much nicer but I'm not sure that code written using pthreads and C++11 threads can
    interoperate.  If not, then we have to make the switch in many libraries at once!)
