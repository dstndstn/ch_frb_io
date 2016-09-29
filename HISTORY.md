- v3: network streams
  - Backwards incompatible: compile flags must now include -pthread
  - Backwards incompatible: re-symlink Makefile.local
  - Describe network streamsx

- v2: implement intensity hdf5 file writing (backwards-compatible update)
  - writing intensity hdf5 files from C++
  - proper unit test of the C++ part
  - new utility ch-plot-intensity-file which makes a quick waterfall plot from an intensity hdf5 file.

- v1: initial version, contains classes for reading intensity hdf5 files and not much else.
