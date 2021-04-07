curDir = $(shell pwd)
CC = g++
CFLAGS  = -Wall -pedantic -O3 -m64 -shared -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 -fPIC

BAMTOOLS_INC = $(curDir)/bamtools/include/bamtools/
BAMTOOLS_LIB = $(curDir)/bamtools/lib/

LREZ_INC = $(curDir)/src/include/
LREZ_LIB = $(curDir)/lib/

LDFLAGS_BAMTOOLS = -lbamtools -L$(BAMTOOLS_LIB)
LDFLAGS_LREZ = -llrez -L$(LREZ_LIB) -lboost_iostreams -lz -lm

MAIN = src/main.o
REVCOMP = src/reverseComplement.o
SOURCE = src/alignmentsRetrieval.o src/barcodesComparison.o src/barcodesExtraction.o src/indexManagementBam.o src/utils.o src/indexManagementFastq.o src/readsRetrieval.o src/gzIndex.o
SUBCOMMANDS = src/subcommands/compare.o src/subcommands/extract.o src/subcommands/help.o src/subcommands/indexBam.o src/subcommands/queryBam.o src/subcommands/indexFastq.o src/subcommands/queryFastq.o

EXEC = bin/LRez
LIB = lib/liblrez.so

all: directories $(LIB) $(EXEC)


directories:
	mkdir -p bin/ lib/

$(LIB): $(SUBCOMMANDS) $(SOURCE) $(REVCOMP)
	$(CC) -fPIC -shared -o $(LIB) $(SUBCOMMANDS) $(SOURCE) $(REVCOMP) $(LDFLAGS_BAMTOOLS) -Wl,-rpath,$(BAMTOOLS_LIB) $(LDFLAGS_LREZ) -Wl,-rpath,$(LREZ_LIB)

$(EXEC): $(MAIN) $(LIB)
	$(CC) -o $(EXEC) $(MAIN) $(LDFLAGS_LREZ) -Wl,-rpath,$(LREZ_LIB)

src/%.o: src/%.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -I$(BAMTOOLS_INC) -I$(LREZ_INC)

src/reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS) -I$(LREZ_INC)


clean:
	rm src/*.o src/subcommands/*.o $(EXEC) $(LIB)


.PHONY: all clean directories
