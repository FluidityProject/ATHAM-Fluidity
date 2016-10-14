SHELL = /bin/sh

#
#  Makefile for PRC code:
#

# Set Makefile's default goal
#
all: libLWM.a

# Set f90 compiler
#
include ../../../Makefile.COMPILE

# Set include path
#

MODDIR =../modules/
CURDIR =./

IFLAGS =$(IGLOBAL) $(MODFLAG)$(CURDIR)$(MODGLUE)$(MODDIR)

# Collect set of option flags included in every compile
#
FLAGS = $(PFLAGS) $(OFLAGS) $(IFLAGS)

# Set flags needed by Makefile.COMMON
# Set the commands to use
#
RM=rm -f 
CP=cp -f 

# Include non-os specifi stuff
#
include Makefile.COMMON

# Rules for compiling f90 files 
#

MOD_MAKE:
	cd $(MODDIR); make

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(F90) -c $(FLAGS) $<

# Generate library
#
libLWM.a: MOD_MAKE $(OBJECTS)
	ar -rs libLWM.a $(ALLOBJECTS)

### DO NOT remove this line - make depends on it ###
