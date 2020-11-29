# Usage:
# make        # compile all binaries
# make clean  # remove ALL binaries and objects

.PHONY	= all clean
RMO			= rm -vf
RMF     = rm -vfr
BINDIR  = ./bin/
MKBINDIR		= mkdir $(BINDIR)

# Compiler
CXX				= g++
CC				= $(CXX)
CXXFLAGS 	= -std=c++11

# Linker
LDFLAGS		:= -L/usr/lib/x86_64-linux-gnu/ -lboost_program_options -pthread

# BlooMine source
SOURCES		= ./CPP_src/BlooMine.cpp
OBJECTS		= $(SOURCES:.cpp=.o)
BINS 			= ./BlooMine

all: $(SOURCES) $(BINS)

$(BINS): $(OBJECTS)
	@echo "Creating object..."
	$(CC) $< -o $@ $(LDFLAGS)

.cpp.o:
	@echo "Checking.."
	$(CC) -c $(CXXFLAGS) $< -o $@

install:
	@echo "Installing..."
	$(MKBINDIR)
	mv $(BINS) $(BINDIR)

clean:
	@echo "Cleaning up..."
	$(RMO) $(OBJECTS)
	$(RMF) ./bin/
