ROOTCONFIG   := root-config
ROOTCINT     := rootcint
ARCH         := $(shell $(ROOTCONFIG) --arch)
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags) -lSpectrum
ROOTGLIBS    := $(shell $(ROOTCONFIG) --libs)

CXX           = g++
CXXFLAGS      = -O3 -Wall -fPIC
LD            = g++
LDFLAGS       = -O3
SOFLAGS       = -shared

INCDIR        =../include

ifeq ($(ARCH),macosx64)
FORTRAN       = gfortran
else
FORTRAN       = g++
endif
FFLAGS        = -fPIC -m64 -pthread
FLFLAGS       = -lgfortran

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS) $(ROOTGLIBS)

CXXFLAGS     += -I$(INCDIR)

NMRANAO		  = NMRana.o
NMRANA        = NMRana

OBJS          = $(NMRANAO)
PROGRAMS      = $(NMRANA) 

all:            $(PROGRAMS)

.SUFFIXES: .cxx .o


$(NMRANA):   $(NMRANAO)
	$(LD) $^ -o $@ $(LDFLAGS)  
	@echo "$@ done."

.SUFFIXES: .cxx

.cxx.o:
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $<

.PHONY: clean

clean:
	@echo "Cleanning everything ... "
	@rm $(PROGRAMS) $(OBJS) $(SLIBS) 
