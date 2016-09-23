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

  - We plan to implement RPC's which operate on the assembled_ringbufs, and support
    operations such as flushes to disk, and retrieval of ring-buffered data over the network.
    The first step here is doing some research to decide which RPC framework to use (maybe
    google RPC?)

  - We don't have thread-safe logging, so diagnostic messages from the various threads
    sometimes interleave each other and are unreadable.  (Implementing thread-safe
    logging is a general todo item for the whole CHIMEFRB backend, not just this repository!)

  - There is a proposal to bitshuffle-compress the packets, which may reduce bandwidth
    by ~20%, but this is currently unimplemented.

  - Eventually it would be nice to do an end-to-end test of the FRB backend, by having
    a "simulator" node generate timestreams containing noise + FRB's, and sending them
    over the network.  However, the current simulation code is single-threaded, and a
    single core is too slow to even generate a stream of Gaussian random numbers at
    full CHIME bandwidth!

  - When the network stream is running, it maintains "event counts" for many types of
    events, such as packet drops, assembler hits/misses etc.  Right now we don't really
    do anything with this information but it's intended to be a starting point for some
    sort of RPC-driven dashboard which can give a visual summary of how the backend is
    performing.

    Related: we probably want to generalize the event counts (currently cumulative) to
    keep track of the event rate for some choice of timescale (say 10 sec).

  - Currently, we have to run `test-network-streams.cpp` at very low throughput (0.1 Gbps)
    to avoid dropping packets.  This means that the unit tests take about an hour to run,
    which isn't really a problem, but is indicative of deeper performance problems?  It
    would be nice to understand where the bottleneck is.

  - It might be possible to further optimize the assembly-language kernels for packet assembly
    and decoding, using streaming writes.  See Chapter 7 ("optimizing cache usage") of the
    intel optimization manual.  We also might be able to improve performance a little using
    aligned loads/stores, but this would impose pointer alignment requirements on callers of
    assembled_chunk::decode().

  - It would be great to switch from pthreads to C++11 threads.  (The C++11 API is
    much nicer but I'm not sure that code written using pthreads and C++11 threads can
    interoperate.  If not, then we have to make the switch in many libraries at once!)

  - Nuisance issue: if a chime_network_stream is constructed from python, then it doesn't
    respond to control-C (not sure if this is a 'ch_frb_io' loose end, or an 'rf_pipelines' 
    loose end).

  - We'll probably want a member function of intensity_network_stream which pins the
    network and assembler threads to specific cores.

  - Open-ended item: there are lots of things that can go wrong in a realtime system,
    such as temporary network failures, and threads running slow so that ring buffers
    overfill.  We need to think carefully about different failure modes and figure out
    how best to handle them.


