#include "subcommands/queryBam.h"
#include "subcommands/help.h"
#include "alignmentsRetrieval.h"
#include "indexManagementBam.h"
#include <getopt.h>

void subcommandQueryBam(int argc, char* argv[]) {
	string bamFile;
	string indexFile;
	string query;
	string list;
	string outputFile;
	BarcodesOffsetsIndex BarcodesOffsetsIndex;
	ofstream out;
	bool outputHeader = false;
	unsigned nbThreads = 1;

	const struct option longopts[] = {
		{"bam",				required_argument,	0, 'b'},
		{"index",			required_argument,	0, 'i'},
		{"query",			required_argument,	0, 'q'},
		{"list",			required_argument,	0, 'l'},
		{"output",			required_argument,	0, 'o'},
		{"threads",			required_argument,	0, 't'},
		{"header",			no_argument,		0, 'H'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "b:i:q:l:o:Ht:", longopts, &index);
	if (iarg == -1) {
		subcommandHelp("query bam");
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'b':
				bamFile = optarg;
				break;
			case 'i':
				indexFile = optarg;
				break;
			case 'q':
				query = optarg;
				break;
			case 'l':
				list = optarg;
				break;
			case 'o':
				outputFile = optarg;
				break;
			case 'H':
				outputHeader = true;
				break;
			case 't':
				nbThreads = static_cast<uint32_t>(stoul(optarg));
				break;
			default:
				subcommandHelp("query bam");
				break;
		}
		iarg = getopt_long(argc, argv, "b:i:q:l:o:Ht:", longopts, &index);
	}

	if (bamFile.empty()) {
		fprintf(stderr, "Please specify a BAM file with option --bam/-b.\n");
		exit(EXIT_FAILURE);
	}
	if (indexFile.empty()) {
		fprintf(stderr, "Please specify an index file with option --index/-i.\n");
		exit(EXIT_FAILURE);
	}
	if (query.empty() and list.empty()) {
		fprintf(stderr, "Please specify a query barcode with option --query/-q or a list of barcodes with option --list/-l.\n");
		exit(EXIT_FAILURE);
	}
	if (!query.empty() and !list.empty()) {
		fprintf(stderr, "Options --query/-q and --list/-l cannot be set at the same time.\n");
		exit(EXIT_FAILURE);
	}

	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	try {
	    BarcodesOffsetsIndex = loadBarcodesOffsetsIndex(indexFile);
		vector<BamAlignment> alignments;
		if (!query.empty()) {
			alignments = retrieveAlignmentsWithBarcode_BamReader(reader, BarcodesOffsetsIndex, query);
		} else {
			alignments = retrieveAlignmentsWithBarcodes(bamFile, BarcodesOffsetsIndex, list, nbThreads);
		}

		if (!outputFile.empty()) {
			out.open(outputFile, ios::out);
			if (!out.is_open()) {
				fprintf(stderr, "Unable to open file %s.", outputFile.c_str());
				exit(EXIT_FAILURE);
			}
			if (outputHeader) {
				out << reader.GetHeaderText() << endl;
			}
			for (BamAlignment al : alignments) {
				out << convertToSam(al, reader.GetReferenceData()) << endl;
			}
			out.close();
		} else {
			if (outputHeader) {
				cout << reader.GetHeaderText() << endl;
			}
			for (BamAlignment al : alignments) {
				cout << convertToSam(al, reader.GetReferenceData()) << endl;
			}
		}

	} catch (exception const& e) {
		cerr << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	reader.Close();
}