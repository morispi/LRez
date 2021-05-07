#include "subcommands/indexFastq.h"
#include "subcommands/help.h"
#include "indexManagementFastq.h"
#include "gzIndex.h"
#include <getopt.h>
#include <chrono>


using namespace std::chrono;

void subcommandIndexFastq(int argc, char* argv[]) {
	string fastqFile;
	string output;
	bool gzipped = false;

	const struct option longopts[] = {
		{"fastq",				required_argument,	0, 'f'},
		{"output",				required_argument,	0, 'o'},
		{"gzip",				no_argument,		0, 'g'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "f:o:g", longopts, &index);
	if (iarg == -1) {
		subcommandHelp("index fastq");
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'f':
				fastqFile = optarg;
				break;
			case 'o':
				output = optarg;
				break;
			case 'g':
				gzipped = true;
				break;
			default:
				subcommandHelp("index fastq");
				break;
		}
		iarg = getopt_long(argc, argv, "f:o:g", longopts, &index);
	}

	if (fastqFile.empty()) {
		fprintf(stderr, "Please specify a fastq file with option --fastq/-f.\n");
		exit(EXIT_FAILURE);
	}
	if (output.empty()) {
		fprintf(stderr, "Please specify an output file with option --output/-o.\n");
		exit(EXIT_FAILURE);
	}

	BarcodesIndex barcodesIndex;
	if (!gzipped) {
		barcodesIndex = indexBarcodesFromFastq(fastqFile);
	} else {
		barcodesIndex = indexBarcodesFromFastqGz(fastqFile);

		// Build and serialize the gzip index
		struct access* index;
		int len = buildGzIndex(fastqFile, SPAN, &index);
		if (len < 0) {
		    switch (len) {
		    case Z_MEM_ERROR:
		        fprintf(stderr, "zran: out of memory\n");
				exit(EXIT_FAILURE);
		    case Z_DATA_ERROR:
		        fprintf(stderr, "zran: compressed data error in %s\n", fastqFile.c_str());
		        exit(EXIT_FAILURE);
		    case Z_ERRNO:
		        fprintf(stderr, "zran: read error on %s\n", fastqFile.c_str());
		        exit(EXIT_FAILURE);
		    default:
		        fprintf(stderr, "zran: error %d while building index\n", len);
		        exit(EXIT_FAILURE);
		    }
		}
		serializeGzIndex(index, fastqFile + "i");
		freeGzIndex(index);
	}
	saveBarcodesIndex(barcodesIndex, output);
}