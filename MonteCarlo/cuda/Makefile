# Compilatori
GCC = gcc
CC = nvc
NVCC = nvcc

# Flags
CFLAGS = -O3 -std=c99 -g -Wall -I../include
NVCCFLAGS = -O3 -g -arch=sm_80 -Xcompiler="-Wall" -I../include
LDFLAGS = -lm -lcudart

# Sorgenti comuni
RANDOM = ranlxs ranlxd gauss
START = start utils
EXTRAS = sweep_ising_2D dH_ising_2D H_ising_2D
CHECK = geometrize_corr_fun geometrize_V geometrize_T

MODULES = $(RANDOM) $(START) $(EXTRAS) $(CHECK)
MDIR = ../modules
VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/extras:$(MDIR)/check

OBJECTS = $(addsuffix .o, $(MODULES))

# Lista degli eseguibili
EXECUTABLES = cuda_Ising_checkerboard cuda_Ising_2D.cu Ising_2d.cu 

# Regola di default
all: $(EXECUTABLES)

# Regola per ogni eseguibile (automaticamente usa .cu se esiste, altrimenti .c)
#%: %.o $(OBJECTS)
#	$(NVCC) $^ -o $@ $(LDFLAGS)

# Compilazione moduli C
%.o: %.c
	$(GCC) $(CFLAGS) -c $< -o $@

# Compilazione moduli CUDA
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Pulizia
clean:
	rm -f *.o $(EXECUTABLES)