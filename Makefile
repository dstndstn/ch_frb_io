# Makefile.local must define the following variables
#   BINDIR       install dir for executables
#   LIBDIR       install dir for C++ libraries
#   INCDIR       install dir for C++ headers
#   CPP          c++ compiler command line (must support c++11)
#   CPP_LFLAGS   optional: any c++ flags which should only be specified when linking
#
# See examples in site/ directory.

include Makefile.local

ifndef CPP
$(error Fatal: Makefile.local must define CPP variable)
endif

ifndef BINDIR
$(error Fatal: Makefile.local must define BINDIR variable)
endif

ifndef INCDIR
$(error Fatal: Makefile.local must define INCDIR variable)
endif

ifndef LIBDIR
$(error Fatal: Makefile.local must define LIBDIR variable)
endif


####################################################################################################


LIBS = -lhdf5

OFILES = assembled_chunk.o \
	assembled_chunk_ringbuf.o \
	avx2_kernels.o \
	hdf5.o \
	intensity_hdf5_file.o \
	intensity_hdf5_ofile.o \
	intensity_network_stream.o \
	intensity_network_ostream.o \
	intensity_packet.o \
	lexical_cast.o \
	udp_packet_list.o \
	udp_packet_ringbuf.o

INCFILES=ch_frb_io.hpp ch_frb_io_internals.hpp
LIBFILES=libch_frb_io.so
INSTALLED_BINARIES=ch-show-intensity-file
INSTALLED_SCRIPTS=ch-plot-intensity-file

TEST_BINARIES = test-intensity-hdf5-file \
	test-misc \
	test-network-streams

all: $(INSTALLED_BINARIES) $(TEST_BINARIES) $(LIBFILES)

install: $(INCFILES) $(LIBFILES) $(INSTALLED_BINARIES)
	mkdir -p $(INCDIR) $(LIBDIR) $(BINDIR)
	cp -f $(INCFILES) $(INCDIR)/
	cp -f $(LIBFILES) $(LIBDIR)/
	cp -f $(INSTALLED_BINARIES) $(INSTALLED_SCRIPTS) $(BINDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(LIBFILES); do rm -f $(LIBDIR)/$$f; done
	for f in $(INSTALLED_BINARIES) $(INSTALLED_SCRIPTS); do rm -f $(BINDIR)/$$f; done

clean:
	rm -f *~ site/*~ *.o *.so $(TEST_BINARIES) $(INSTALLED_BINARIES)

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

libch_frb_io.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -o $@ -shared $^ $(LIBS)

ch-show-intensity-file: ch-show-intensity-file.cpp $(INCFILES) libch_frb_io.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lch_frb_io

test-intensity-hdf5-file: test-intensity-hdf5-file.cpp $(INCFILES) libch_frb_io.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lch_frb_io

test-misc: test-misc.cpp $(INCFILES) libch_frb_io.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lch_frb_io

test-network-streams: test-network-streams.cpp $(INCFILES) libch_frb_io.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lch_frb_io

ringbuf: ringbuf.cpp $(OFILES)
	$(CPP) $(CPP_LFLAGS) -o $@ $^ -lhdf5
