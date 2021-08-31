#include "subcommands/compare.h"
#include "subcommands/extract.h"
#include "subcommands/help.h"
#include "subcommands/indexBam.h"
#include "subcommands/queryBam.h"
#include "subcommands/indexFastq.h"
#include "subcommands/queryFastq.h"
#include "subcommands/stats.h"
#include <set>
#include <string.h>

using namespace std;


int main(int argc, char* argv[]) {
	set<string> subcommands {"index", "query", "extract", "compare", "stats"};
	if (argc < 2 or !subcommands.count(argv[1])) {
		subcommandHelp("global");
	} else if (!strcmp(argv[1], "compare")) {
		subcommandCompare(argc, argv);
	} else if (!strcmp(argv[1], "extract")) {
		subcommandExtract(argc, argv);
	} else if (!strcmp(argv[1], "stats")) {
		subcommandStats(argc, argv);
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
