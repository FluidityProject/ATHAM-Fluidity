SHELL = /bin/sh

#
#  Makefile for PRC code:
#

# Set Makefile's default goal
#
all: make_all_objects

# Set f90 compiler
#
include ../../../Makefile.COMPILE

# Set include path
#

all: make_all_objects

IFLAGS =$(IGLOBAL) 

# Collect set of option flags included in every compile
#
FLAGS = $(PFLAGS) $(OFLAGS) $(IFLAGS)

# Set flags needed by Makefile.COMMON
# Set the commands to use
#
RM=rm -f 
CP=cp -f 

# Include non-os specific stuff
#
include Makefile.COMMON

# Rules for compiling f90 files 
#
.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(F90) -c $(FLAGS) $<

# Generate library
#
make_all_objects: $(MOD_OBJECTS)
#	ar -rs libLWM.a $(ALLOBJECTS)

### DO NOT remove this line - make depends on it ###
