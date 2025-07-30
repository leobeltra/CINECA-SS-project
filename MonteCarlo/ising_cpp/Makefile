# Compiler settings
CXX = g++
NVCC = nvcc
CXXFLAGS = -Wall -I../include -I/usr/local/cuda-12.9/include
NVCCFLAGS = -arch=sm_89 -I../include

# Output executable
TARGET_BASE = naive tiles2d

# Source files
CPP_SOURCES = main.cpp cpu.cpp $(wildcard utils/*.cpp)

# Default target
TARGET = $(addsuffix .exe,$(TARGET_BASE))
all: $(TARGET)

# Link rule
# Define function to get objects for a target
target_objects = $(patsubst %.cpp,%.o,$(wildcard $1/*.cpp)) $(patsubst %.cu,%.o,$(wildcard $1/*.cu))

# Generate target objects for each base target
all_targets_objects := $(foreach dir,$(TARGET_BASE),$(call target_objects,$(dir)))

# Link rule
%.exe: $(CPP_SOURCES:.cpp=.o) $(all_targets_objects)
	@echo "Linking $@ with sources from $* directory"
	@echo "Target objects: $(call target_objects,$*)"
	$(NVCC) $(NVCCFLAGS) $(CPP_SOURCES:.cpp=.o) $(call target_objects,$*) -lm -lcudart -o $@

# Compile C++ sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA sources
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f ./*.o $(TARGET) $(all_targets_objects)

.PHONY: all clean