#include "subcommands/queryFastq.h"
#include "subcommands/help.h"
#include "readsRetrieval.h"
#include "indexManagementFastq.h"
#include <getopt.h>

void subcommandQueryFastq(int argc, char* argv[]) {
	string fastqFile;
	string indexFile;
	string query;
	string list;
	string outputFile;
	BarcodesIndex BarcodesIndex;
	BamReader reader;
	ofstream out;
	bool gzipped = false;

	const struct option longopts[] = {
		{"fastq",			required_argument,	0, 'f'},
		{"index",			required_argument,	0, 'i'},
		{"query",			required_argument,	0, 'q'},
		{"list",			required_argument,	0, 'l'},
		{"output",			required_argument,	0, 'o'},
		{"gzipped",			no_argument,		0, 'g'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "f:i:q:l:o:g", longopts, &index);
	if (iarg == -1) {
		subcommandHelp("query fastq");
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'f':
				fastqFile = optarg;
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
			case 'g':
				gzipped = true;
				break;
			default:
				subcommandHelp("query bam");
				break;
		}
		iarg = getopt_long(argc, argv, "f:i:q:l:o:g", longopts, &index);
	}

	if (fastqFile.empty()) {
		fprintf(stderr, "Please specify a fastq file with option --fastq/-f.\n");
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

	BarcodesIndex = loadBarcodesIndex(indexFile);
	vector<string> reads;
	if (!query.empty()) {
		if (!gzipped) {
			reads = retrieveReadsWithBarcode(fastqFile, BarcodesIndex, query);
		} else {
			reads = retrieveReadsWithBarcode_Gzip(fastqFile, BarcodesIndex, query);
		}
	} else {
		if (!gzipped) {
			reads = retrieveReadsWithBarcodes(fastqFile, BarcodesIndex, list);
		} else {
			reads = retrieveReadsWithBarcodes_Gzip(fastqFile, BarcodesIndex, list);
		}
	}

	if (!outputFile.empty()) {
		out.open(outputFile, ios::out);
		if (!out.is_open()) {
			fprintf(stderr, "Unable to open file %s.", outputFile.c_str());
			exit(EXIT_FAILURE);
		}
		for (string r : reads) {
			out << r << endl;
		}
		out.close();
	} else {
		for (string r : reads) {
			cout << r << endl;
		}
	}
}