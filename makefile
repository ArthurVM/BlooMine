# Usage:
# make        # compile all binaries
# make clean  # remove ALL binaries and objects

.PHONY	= all clean
RMO			= rm -vf
RMF     = rm -vfr
MKDIR		= mkdir ./bin/

# Compiler
CXX				= g++
CC				= $(CXX)
CXXFLAGS 	= -std=c++11

# Linker
LDFLAGS		:= -L/usr/lib/x86_64-linux-gnu/ -lboost_program_options -pthread

# BlooMine source
SOURCES		= ./CPP_src/BlooMine.cpp
OBJECTS		= $(SOURCES:.cpp=.o)
BINS 			= ./bin/BlooMine

all: $(SOURCES) $(BINS)

$(BINS): $(OBJECTS)
	@echo "Making bin dir."
	$(MKDIR)
	@echo "Creating object..."
	$(CC) $< -o $@ $(LDFLAGS)

.cpp.o:
	@echo "Checking.."
	$(CC) -c $(CXXFLAGS) $< -o $@

clean:
	@echo "Cleaning up..."
	$(RMO) $(OBJECTS)
	$(RMF) ./bin/
