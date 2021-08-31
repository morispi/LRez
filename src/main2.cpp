#include "subcommands/compare.h"
#include "subcommands/extract.h"
#include "subcommands/help.h"
#include "subcommands/indexBam.h"
#include "subcommands/queryBam.h"
#include "subcommands/indexFastq.h"
#include "subcommands/queryFastq.h"
#include "barcodesLoading.h"
#include <set>
#include <string.h>

using namespace std;

// Path were the barcodes files are stored
string path;


int main(int argc, char* argv[]) {
	// Retrieve path to the barcodes folder
	string tmp(argv[0]);
	path = tmp.substr(0, tmp.length() - 4) + "..share/--/";

	set<string> subcommands {"index", "query", "extract", "compare"};
	if (argc < 2 or !subcommands.count(argv[1])) {
		subcommandHelp("global");
	} else if (!strcmp(argv[1], "compare")) {
		subcommandCompare(argc, argv);
	} else if (!strcmp(argv[1], "extract")) {
		subcommandExtract(argc, argv);
	} else if (!strcmp(argv[1], "index")) {
		if (argc < 3) {
			subcommandHelp("global");
		}
		if (!strcmp(argv[2], "bam")) {
			subcommandIndexBam(argc, argv);
		} else if (!strcmp(argv[2], "fastq")) {
			subcommandIndexFastq(argc, argv);
		} else {
			subcommandHelp("global");
		}
	} else if (!strcmp(argv[1], "query")) {
		if (argc < 3) {
			subcommandHelp("global");
		}
		if (!strcmp(argv[2], "bam")) {
			subcommandQueryBam(argc, argv);
		} else if (!strcmp(argv[2], "fastq")) {
			subcommandQueryFastq(argc, argv);
		} else {
			subcommandHelp("global");
		}
	} else {
		fprintf(stderr, "An unexpected error happened. Please try running again.\n");
		exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;
}
