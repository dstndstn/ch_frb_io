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


INCFILES=ch_frb_io.hpp ch_frb_io_internals.hpp
LIBFILES=libch_frb_io.so
INSTALLED_BINARIES=ch-show-intensity-file

LIBS=-lhdf5
OFILES=hdf5.o intensity_hdf5_file.o

all: $(INSTALLED_BINARIES) $(LIBFILES)

install: $(INCFILES) $(LIBFILES) $(INSTALLED_BINARIES)
	mkdir -p $(INCDIR) $(LIBDIR) $(BINDIR)
	cp -f $(INCFILES) $(INCDIR)/
	cp -f $(LIBFILES) $(LIBDIR)/
	cp -f $(INSTALLED_BINARIES) $(BINDIR)/

uninstall:
	for f in $(INCFILES); do rm -f $(INCDIR)/$$f; done
	for f in $(LIBFILES); do rm -f $(LIBDIR)/$$f; done
	for f in $(BINFILES); do rm -f $(BINDIR)/$$f; done

clean:
	rm -f *~ site/*~ *.o *.so $(INSTALLED_BINARIES)

%.o: %.cpp $(INCFILES)
	$(CPP) -c -o $@ $<

libch_frb_io.so: $(OFILES)
	$(CPP) $(CPP_LFLAGS) -o $@ -shared $^ $(LIBS)

ch-show-intensity-file: ch-show-intensity-file.cpp $(INCFILES) libch_frb_io.so
	$(CPP) $(CPP_LFLAGS) -o $@ $< -lch_frb_io
