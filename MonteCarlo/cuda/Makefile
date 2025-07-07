# CUDA Makefile using nvc for .cu

# Compilatori
GCC = gcc         
CC = nvc
NVCC = nvcc

# Flags per file .c
CFLAGS = -O2 -std=c99 -g -Wall -I../include

# Flags per file .cu
NVCCFLAGS = -O2 -g -arch=sm_80 -Xcompiler "-Wall" -I../include
LDFLAGS = -lm -lcudart

RANDOM = ranlxs ranlxd gauss

START = start utils

EXTRAS = dazione_an azione_anharmonic sweep_an_creutz dazione_an_creutz azione_anharmonic_creutz dazione_2d dazione_2d_an sweep_2d d_azione_klauder dH_ising H_Ising sweep_ising sweep_ising_2D dH_ising_2D H_ising_2D

MONTECARLO = azione dazione sweep

CHECK = geometrize_corr_fun geometrize_V geometrize_T

JACKKNIFE = 

SCALAR_FIELD_THEORY = sweep_scalar_field scalar_field_action

BINNING = binned_correlation correlation

MODULES = $(RANDOM) $(START) $(MONTECARLO) $(EXTRAS) $(CHECK) $(JACKKNIFE) $(SCALAR_FIELD_THEORY)

# search path for modules

MDIR = ../modules

VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/extras:$(MDIR)/montecarlo:$(MDIR)/check:$(MDIR)/jackknife:$(MDIR)/scalar_field_theory

OBJECTS = $(addsuffix .o, $(MODULES))

TARGET = cuda_Ising_2D

all: $(TARGET)

# Compilazione moduli C con GCC
%.o: %.c
	$(GCC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(TARGET): $(TARGET).o $(OBJECTS)
	$(NVCC) $(NVCCFLAGS) $(TARGET).o $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm -f *.o *.mod $(TARGET)
