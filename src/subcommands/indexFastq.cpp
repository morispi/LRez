#include "subcommands/indexFastq.h"
#include "subcommands/help.h"
#include "indexManagementFastq.h"
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
		{"userx",				no_argument,		0, 'u'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "f:o:gu", longopts, &index);
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
			case 'u':
				CONSIDER_RX = true;
				break;
			default:
				subcommandHelp("index fastq");
				break;
		}
		iarg = getopt_long(argc, argv, "f:o:gu", longopts, &index);
	}

	if (fastqFile.empty()) {
		fprintf(stderr, "Please specify a fastq file with option --fastq/-f.\n");
		exit(EXIT_FAILURE);
	}
	if (output.empty()) {
		fprintf(stderr, "Please specify an output file with option --output/-o.\n");
		exit(EXIT_FAILURE);
	}

	cerr << "parse argument" << endl;
	BarcodesIndex barcodesIndex;
	barcodesIndex = indexBarcodesFromFastq(fastqFile);
	saveBarcodesIndex(barcodesIndex, output);
}