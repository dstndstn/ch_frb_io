- v3: large update containing an intial implementation of chimefrb networking code.

  - Backward incompatible: the C++ source files have been moved from cpp/ to the
    top-level directory, for consistency with other chimefrb repos.  This is a 
    "backward incompatible" change only in the sense that your Makefile.local
    should now be in the toplevel directory, not in cpp/.

  - Backwards incompatible: compiler flags must now include -pthread

  - Networking code has been written.  The most important classes are:

      - intensity_network_stream, which receives a packet stream over the 
        network, decodes it and assembles a sequence of regular arrays for 
        each beam.

      - intensity_network_ostream, which receives a sequence of regular arrays,
        encodes and packetizes the data, and sends it over the network.

- v2: implement intensity hdf5 file writing (backwards-compatible update)
  - writing intensity hdf5 files from C++
  - proper unit test of the C++ part
  - new utility ch-plot-intensity-file which makes a quick waterfall plot from an intensity hdf5 file.

- v1: initial version, contains classes for reading intensity hdf5 files and not much else.
