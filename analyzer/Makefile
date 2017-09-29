# -*- Makefile -*- (for emacs)
# src/Makefile-UserModules.  Generated from Makefile-UserModules.in by configure.
#
# Herwig 7.1.1
#
# This Makefile is intended for compiling Herwig plugins
# You can find example plugins in the Contrib directory
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

#
# C++ flags
# 
CXX = g++ -std=c++11

CXXFLAGS = -g -W -Wall -O2

# add additional package header directories here
CPPFLAGS = -I/mnt/t3nfs01/data01/shome/theofil/GEN/./include -I/mnt/t3nfs01/data01/shome/theofil/GEN/./include -I/mnt/t3nfs01/data01/shome/theofil/GEN/./include

ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS  := $(shell root-config --glibs)
ROOT       = $(ROOTCFLAGS) $(ROOTGLIBS)


# add additional libraries here
LDFLAGS = 

ALLCCFILES=$(shell echo *.cc)
all : $(ALLCCFILES:.cc=.so)

#%.so : %.cc %.h 
#	$(CXX) -shared -fPIC $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< -o $@  
%.so : %.cc %.h 
	$(CXX) -shared -fPIC $(CPPFLAGS) $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) $< -o $@  

clean:
	rm -f $(ALLCCFILES:.cc=.so)
