#include "subcommands/extract.h"
#include "subcommands/help.h"
#include "barcodesExtraction.h"
#include <getopt.h>

using namespace std;

void subcommandExtract(int argc, char* argv[]) {
	string bamFile;
	string region;
	bool extractAll = false;
	bool duplicates = false;
	string outputFile;
	BamReader reader;
	ofstream out;
	BamRegion r;
	unsigned nbThreads = 1;

	const struct option longopts[] = {
		{"bam",			required_argument,	0, 'b'},
		{"region",		required_argument,	0, 'r'},
		{"all",			no_argument,		0, 'a'},
		{"duplicates",	no_argument,		0, 'd'},
		{"output",		required_argument,	0, 'o'},
		{"threads",		required_argument,	0, 't'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "b:r:o:adt:", longopts, &index);
	if (iarg == -1) {
		subcommandHelp("extract");
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'b':
				bamFile = optarg;
				break;
			case 'r':
				region = optarg;
				break;
			case 'a':
				extractAll = true;
				break;
			case 'o':
				outputFile = optarg;
				break;
			case 'd':
				duplicates = true;
				break;
			case 't':
				nbThreads = static_cast<uint32_t>(stoul(optarg));
				break;
			default:
				subcommandHelp("extract");
				break;
		}
		iarg = getopt_long(argc, argv, "b:r:o:adt:", longopts, &index);
	}

	if (bamFile.empty()) {
		fprintf(stderr, "Please specify a BAM file with option --bam/-b.\n");
		exit(EXIT_FAILURE);
	}
	if (!extractAll and region.empty()) {
		fprintf(stderr, "Please specify a region of interest with option --region/-r, or use option --all/-a to extract all barcodes.\n");
		exit(EXIT_FAILURE);
	}
	if (extractAll and !region.empty()) {
		fprintf(stderr, "Options --region/-r and --all/-a cannot be used at the same time.\n");
		exit(EXIT_FAILURE);
	}

	try {
		if (!duplicates) {
			robin_hood::unordered_set<string> barcodes;

			if (extractAll) {
				barcodes = extractBarcodesSeqsFromBAM(bamFile, nbThreads);
			} else {
				barcodes = extractBarcodesSeqsFromRegion(bamFile, region);
			}

			if (!outputFile.empty()) {
				out.open(outputFile, ios::out);
				if (!out.is_open()) {
					fprintf(stderr, "Unable to open file %s.", outputFile.c_str());
					exit(EXIT_FAILURE);
				}
				for (string b : barcodes) {
					out << b << endl;
				}
				out.close();
			} else {
				for (string b : barcodes) {
					cout << b << endl;
				}
			}		
		} else {
			vector<string> barcodes;

			if (extractAll) {
				barcodes = extractBarcodesSeqsFromBAMWithDuplicates(bamFile, nbThreads);
			} else {
				barcodes = extractBarcodesSeqsFromRegionWithDuplicates(bamFile, region);
			}

			if (!outputFile.empty()) {
				out.open(outputFile, ios::out);
				if (!out.is_open()) {
					fprintf(stderr, "Unable to open file %s.", outputFile.c_str());
					exit(EXIT_FAILURE);
				}
				for (string b : barcodes) {
					out << b << endl;
				}
				out.close();
			} else {
				for (string b : barcodes) {
					cout << b << endl;
				}
			}
		}
	} catch (exception const& e) {
		cerr << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	
}