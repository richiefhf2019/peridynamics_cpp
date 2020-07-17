# options available in this makefile:
# 1. option to toggle on/off OpenMP (default on)
# 2. option to toggle on/off GCC quadmath (default off)
# 3. option to switch between optimization/debug mode (default optimization)

# C++ compiler
CXX = g++

# optimization
OPTIMIZE = -O3 -DNDEBUG

# debugging options
DEBUG = -Wall -g -pg

# 0. C++11
C11 = -std=c++0x

# 1. OpenMP (macro OPENMP is defined in source code)
OPENMP = -DOPENMP -fopenmp
#OPENMP =

# 2. GCC quadmath
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# 3. Optimize or debug
#CXXFLAGS = $(OPTIMIZE) $(C11) $(OPENMP) $(QUADMATH)
CXXFLAGS = $(DEBUG) $(QUADMATH)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = periDynamics

.PHONY: all tar clean

# template classes contain implementation and need to recompile upon change
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

sinclude $(SOURCES:.cpp=.d)

%.d: %.cpp
	$(CXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

tar:
	tar -cvf $(EXECUTABLE).tar *.h *.cpp makefile* readme
clean:
	-rm -f *.o *.d  *~ *.tar $(EXECUTABLE)

# sinclude is always resolved even if make tar/clean
