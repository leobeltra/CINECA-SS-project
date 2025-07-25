# Compiler settings
CXX = g++
NVCC = nvcc
CXXFLAGS = -Wall -I../include -I/usr/local/cuda-12.9/include
NVCCFLAGS = -arch=sm_89 -I../include

# Output executable
TARGET = Ising.exe

# Source files
CPP_SOURCES = cpu.cpp main.cpp
CU_SOURCES = gpu.cu
OBJECTS = $(CPP_SOURCES:.cpp=.o) $(CU_SOURCES:.cu=.o)

# Default target
all: $(TARGET)

# Link rule
$(TARGET): $(OBJECTS)
	$(NVCC) $(NVCCFLAGS) $^ -lm -lcudart -o $@

# Compile C++ sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA sources
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f ./*.o $(TARGET)

.PHONY: all clean