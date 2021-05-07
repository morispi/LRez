#include "subcommands/indexBam.h"
#include "subcommands/help.h"
#include "indexManagementBam.h"
#include <getopt.h>
#include <chrono>


using namespace std::chrono;

void subcommandIndexBam(int argc, char* argv[]) {
	string bamFile;
	string output;
	bool indexOffsets = false;
	bool indexPositions = false;
	bool primary = false;
	unsigned quality = 0;

	const struct option longopts[] = {
		{"bam",				required_argument,	0, 'b'},
		{"output",			required_argument,	0, 'o'},
		{"offsets",			no_argument,		0, 'f'},
		{"positions",		no_argument,		0, 'p'},
		{"primary",			no_argument,		0, 'r'},
		{"quality",			required_argument,	0, 'q'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "b:o:fprq:", longopts, &index);
	if (iarg == -1) {
		subcommandHelp("index bam");
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'b':
				bamFile = optarg;
				break;
			case 'o':
				output = optarg;
				break;
			case 'f':
				indexOffsets = true;
				break;
			case 'p':
				indexPositions = true;
				break;
			case 'r':
				primary = true;
				break;
			case 'q':
				quality = static_cast<uint32_t>(stoul(optarg));
				break;
			default:
				subcommandHelp("index bam");
				break;
		}
		iarg = getopt_long(argc, argv, "b:o:fprq:", longopts, &index);
	}

	if (bamFile.empty()) {
		fprintf(stderr, "Please specify a BAM file with option --bam/-b.\n");
		exit(EXIT_FAILURE);
	}
	if (output.empty()) {
		fprintf(stderr, "Please specify an output file with option --output/-o.\n");
		exit(EXIT_FAILURE);
	}
	
	if (!indexOffsets and !indexPositions) {
		fprintf(stderr, "Please specify whether you want to index the offsets of the barcodes in the BAM file (--offsets/-f) or the (chromosome, begPosition) pairs of the barcodes (--positions/-p).\n");
		exit(EXIT_FAILURE);
	}

	if (indexOffsets and indexPositions) {
		fprintf(stderr, "Options --offsets/-f and --positions/-p cannot be used at the same time.\n");
		exit(EXIT_FAILURE);
	}

	if (indexOffsets) {
		BarcodesOffsetsIndex barcodesOffsetsIndex;
		barcodesOffsetsIndex = indexBarcodesOffsetsFromBam(bamFile, primary, quality);
		saveBarcodesOffsetsIndex(barcodesOffsetsIndex, output);
	} else if (indexPositions) {
		BarcodesPositionsIndex barcodesPositionsIndex;
		barcodesPositionsIndex = indexBarcodesPositionsFromBam(bamFile, primary, quality);
		saveBarcodesPositionsIndex(barcodesPositionsIndex, output);
	} else {
		fprintf(stderr, "An unexpected error has occurred. Please try running again.\n");
		exit(EXIT_FAILURE);
	}
}