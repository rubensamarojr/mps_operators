# https://latedev.wordpress.com/2014/11/08/generic-makefiles-with-gcc-and-gnu-make/ #
# Executable #
BIN_NAME := myprogram.exe

# Compiler #
CXX := g++
LINKER := g++

# Flags #
INCDIRS := -I. -Iinclude/
CXXFLAGS := -std=c++11 -Wall -Wextra

# Find all source files in the source directory src #
SRCFILES := $(wildcard src/*.cpp) 
OBJFILES := $(patsubst %.cpp,%.o,$(SRCFILES))
DEPFILES := $(patsubst %.cpp,%.d,$(SRCFILES))

$(BIN_NAME): $(OBJFILES)
	$(LINKER) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c $< -o $@

%.d: %.cpp
	$(CXX) $(INCDIRS) -MM $< > $@

-include $(DEPFILES)


