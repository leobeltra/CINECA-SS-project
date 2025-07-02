# CUDA Makefile using nvc for .cu

# Compilatori
GCC = gcc         
CC = nvc
NVCC = nvcc

# Flags per file .c
CFLAGS = -O3 -std=c99 -g -Wall -I../include

# Flags per file .cu
NVCCFLAGS = -O3 -g -arch=sm_80 -Xcompiler="-fopenmp, -Wall" -I../include
LDFLAGS = -lm -lcudart

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

OBJECTS = $(addsuffix .o, $(MODULES))

TARGET = cuda_Ising_checkerboard

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
