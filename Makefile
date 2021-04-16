ifeq ($(SHLIB_EXT),)
	OS_NAME := $(shell uname -s | tr A-Z a-z)
	ifeq ($(OS_NAME),linux)
		SHLIB_EXT=.so
	else
		SHLIB_EXT=.dylib
	endif
endif

curDir = $(shell pwd)

PREFIX ?= /usr/local
BINDIR = $(PREFIX)/bin
LIBDIR = $(PREFIX)/lib

CXX ?= g++
CXXFLAGS += -Wall -pedantic -O3 -m64 -std=c++11 -fPIC
LDFLAGS += -L$(curDir)/lib -Wl,-rpath,$(LIBDIR) -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

BAMTOOLS_LIB_PREFIX = lrez_
BAMTOOLS_LIB = lib/lib$(BAMTOOLS_LIB_PREFIX)bamtools$(SHLIB_EXT)

BAMTOOLS_INC = ./include/bamtools/
LREZ_INC = ./src/include/

LIBS_LREZ = -llrez
LIBS_BOOST_LZ_LM = -lboost_iostreams -lz -lm -lc
LIBS_BAMTOOLS = -l$(BAMTOOLS_LIB_PREFIX)bamtools

MAIN = src/main.o
REVCOMP = src/reverseComplement.o
SOURCE = src/alignmentsRetrieval.o src/barcodesComparison.o src/barcodesExtraction.o src/indexManagementBam.o src/utils.o src/indexManagementFastq.o src/readsRetrieval.o src/gzIndex.o
SUBCOMMANDS = src/subcommands/compare.o src/subcommands/extract.o src/subcommands/help.o src/subcommands/indexBam.o src/subcommands/queryBam.o src/subcommands/indexFastq.o src/subcommands/queryFastq.o

EXEC = bin/LRez
LIB = lib/liblrez${SHLIB_EXT}

all: $(LIB) $(EXEC)

directories:
	mkdir -p bin/ lib/

$(BAMTOOLS_LIB):
	mkdir -p bamtools/build
	cd bamtools/build && \
		cmake \
			-DBUILD_SHARED_LIBS=ON \
			-DCMAKE_SHARED_LIBRARY_PREFIX_CXX=lib$(BAMTOOLS_LIB_PREFIX) \
			-DCMAKE_INSTALL_LIBDIR=lib \
			-DCMAKE_INSTALL_PREFIX=$(curDir) \
			.. && \
		$(MAKE) && \
		$(MAKE) install

$(LIB): $(SUBCOMMANDS) $(SOURCE) $(REVCOMP) $(BAMTOOLS_LIB) directories
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -o $(LIB) $(SUBCOMMANDS) $(SOURCE) $(REVCOMP) $(LIBS_BAMTOOLS) $(LIBS_BOOST_LZ_LM)

$(EXEC): $(MAIN) $(LIB) directories
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $(EXEC) $(MAIN) $(LIBS_LREZ) $(LIBS_BAMTOOLS) $(LIBS_BOOST_LZ_LM)

src/%.o: src/%.cpp $(BAMTOOLS_LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(BAMTOOLS_INC) -I$(LREZ_INC) -o $@ -c $<

src/reverseComplement.o: src/reverseComplement.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(LREZ_INC) -o src/reverseComplement.o -c src/reverseComplement.cpp


install: $(EXEC) $(LIB)
	mkdir -p $(DESTDIR)$(BINDIR)
	mkdir -p $(DESTDIR)$(LIBDIR)
	cp -a $(EXEC) $(DESTDIR)$(BINDIR)/
	cp -a $(LIB)* $(DESTDIR)$(LIBDIR)/
	cp -a $(BAMTOOLS_LIB)* $(DESTDIR)$(LIBDIR)/

clean:
	rm src/*.o src/subcommands/*.o $(EXEC) $(LIB)


.PHONY: all clean directories
