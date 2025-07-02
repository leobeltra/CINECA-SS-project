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

MAIN = prova_azione prova_dazione check_azione check_corr_fun

RANDOM = ranlxs ranlxd gauss

START = start utils

EXTRAS = dazione_an azione_anharmonic sweep_an_creutz dazione_an_creutz azione_anharmonic_creutz dazione_2d dazione_2d_an sweep_2d d_azione_bf azione_bf d_azione_klauder 

MONTECARLO = azione dazione sweep 

CHECK = geometrize_T geometrize_V geometrize_corr_fun

MODULES = $(RANDOM) $(START) $(MONTECARLO) $(EXTRAS) $(CHECK)


# search path for modules

MDIR = ../../modules

VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/extras:$(MDIR)/montecarlo:$(MDIR)/check


# additional include directories

INCPATH = ../../include


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
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


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
