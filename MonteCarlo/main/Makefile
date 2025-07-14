################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all


# main programs and required modules 

MAIN = corr_fun osservabili osservabili_armonic termalization binning_deco ground_state autocorrelation corr_fun_anarmonica ground_state_an_creutz GS_harm_2d GS_dwell corr_fun_Hermite2 autocorrelation_C_binned Ising_1D Ising_2D Spin_Glass scalar_field Ising_2D_serial

RANDOM = ranlxs ranlxd gauss

START = start utils

EXTRAS = sweep_ising_2D dH_ising_2D H_ising_2D

MONTECARLO = 

CHECK = geometrize_corr_fun geometrize_V geometrize_T

JACKKNIFE = 

SCALAR_FIELD_THEORY = 

BINNING = 

MODULES = $(RANDOM) $(START) $(MONTECARLO) $(EXTRAS) $(CHECK) $(JACKKNIFE) $(SCALAR_FIELD_THEORY)


# search path for modules

MDIR = ../modules

VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/extras:$(MDIR)/montecarlo:$(MDIR)/check:$(MDIR)/jackknife:$(MDIR)/scalar_field_theory

# additional include directories

INCPATH = ../include


# additional libraries to be included 
 
LIBS = m

LIBPATH = 


# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -std=c99 -pedantic -fstrict-aliasing \
         -Wall -Wno-long-long -O 
 

############################## do not change ###################################

SHELL=/bin/bash

CC=gcc

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) $< -c -g $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory 

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
