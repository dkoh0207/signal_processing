#
# This is an example GNUmakefile for my packages
#
NAME = sigproc_multithreading

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# include options for this package
INCFLAGS  = -I. -I/sdf/group/neutrino/koh0207/DSP/TBB/tbb/include -I/sdf/group/neutrino/koh0207/DSP/TBB/pstl/stdlib -I/sdf/group/neutrino/koh0207/DSP/TBB/pstl/include -I/sdf/group/neutrino/koh0207/DSP/eigen -I/sdf/group/neutrino/koh0207/DSP/fftw3/include
LDFLAGS   = -L/sdf/group/neutrino/koh0207/DSP/TBB/tbb/lib -L/sdf/group/neutrino/koh0207/DSP/TBB/pstl/stdlib -L/sdf/group/neutrino/koh0207/DSP/fftw3/lib -lfftw3f

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(MYSW_DIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(MYSW_DIR)/Makefile/GNUmakefile.CORE
