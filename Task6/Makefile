CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17 `root-config --cflags`
LDFLAGS = `root-config --libs` -lMinuit2

TARGETS = D_rho levy_fit
SRCS = D_rho.cpp levy_fit.cpp

all: $(TARGETS)

D_rho: D_rho.cpp
	$(CXX) $(CXXFLAGS) -o D_rho D_rho.cpp $(LDFLAGS)

levy_fit: levy_fit.cpp
	$(CXX) $(CXXFLAGS) -o levy_fit levy_fit.cpp $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o
