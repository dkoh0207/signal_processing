#
# This is an example GNUmakefile for my packages
#
NAME = PyUtil

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell python3-config --includes)
INCFLAGS += -I$(shell python3 -c "import numpy; print(numpy.get_include())")

# linker flags
LDFLAGS += -L$(shell python3-config --prefix)/lib/
LDFLAGS += $(shell ${LARCV_PYTHON_CONFIG} --ldflags)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(MYSW_DIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(MYSW_DIR)/Makefile/GNUmakefile.CORE
