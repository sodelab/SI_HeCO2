# ----------------------------------------------- #
# Makefile for CO2He PES (macOS arm64)
# Produces: libCO2He.dylib
# ----------------------------------------------- #

2BRG = CO2He
2B   = mbCO2RG

CXX      = clang++
OPT      = -O2
CXXFLAGS = -fPIC -std=c++11 -Wall -Wextra
LDFLAGS  = -dynamiclib
# Optional: set the install name so loaders can find it relative to @rpath
SONAME   = -Wl,-install_name,@rpath/lib$(2BRG).dylib

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

all: lib$(2BRG).dylib

lib$(2BRG).dylib: $(OBJ)
	$(CXX) $(OPT) $(LDFLAGS) $(SONAME) -o $@ $(OBJ)

clean:
	rm -f *.o *.so *.dylib

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OPT) -c $< -o $@
