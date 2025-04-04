# Makefile to compile, run, and plot Ising model cumulant ratio

CXX = g++
CXXFLAGS = -O3 -fopenmp -march=native -std=c++11
TARGET = ising
SRC = ising_sim.cc
OUTDIR = output
DATAFILE = $(OUTDIR)/ising_sim.txt
PLOT = $(OUTDIR)/ising_sim.png

.PHONY: all clean

all: $(PLOT)

# Generate plot from simulation data
$(PLOT): $(DATAFILE)
	python ising_sim.py

# Run simulation and save output
$(DATAFILE): $(TARGET) | $(OUTDIR)
	./$(TARGET) > $(DATAFILE)

# Compile C++ code
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Create output directory if missing
$(OUTDIR):
	mkdir -p $(OUTDIR)

# Clean generated files
clean:
	rm -f $(TARGET)
	rm -rf $(OUTDIR)
