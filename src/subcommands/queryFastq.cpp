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
	string collectionOfLists;
	string outputFile;
	BarcodesIndex BarcodesIndex;
	BamReader reader;
	ifstream in;
	ofstream out;
	bool gzipped = false;
	unsigned nbThreads = 1;

	const struct option longopts[] = {
		{"fastq",				required_argument,	0, 'f'},
		{"index",				required_argument,	0, 'i'},
		{"query",				required_argument,	0, 'q'},
		{"list",				required_argument,	0, 'l'},
		{"collectionOfLists",	required_argument,	0, 'c'},
		{"output",				required_argument,	0, 'o'},
		{"threads",				required_argument,	0, 't'},
		{"gzipped",				no_argument,		0, 'g'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "f:i:q:l:c:o:gt:", longopts, &index);
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
			case 'c':
				collectionOfLists = optarg;
				break;
			case 'o':
				outputFile = optarg;
				break;
			case 'g':
				gzipped = true;
				break;
			case 't':
				nbThreads = static_cast<uint32_t>(stoul(optarg));
				break;
			default:
				subcommandHelp("query bam");
				break;
		}
		iarg = getopt_long(argc, argv, "f:i:q:l:o:gt:", longopts, &index);
	}

	if (fastqFile.empty()) {
		fprintf(stderr, "Please specify a fastq file with option --fastq/-f.\n");
		exit(EXIT_FAILURE);
	}
	if (indexFile.empty()) {
		fprintf(stderr, "Please specify an index file with option --index/-i.\n");
		exit(EXIT_FAILURE);
	}
	if (query.empty() and list.empty() and collectionOfLists.empty()) {
		fprintf(stderr, "Please specify a query barcode with option --query/-q or a list of barcodes with option --list/-l or a collection of lists of barcodes with option --collectionOfLists/-c.\n");
		exit(EXIT_FAILURE);
	}
	if (!query.empty() and !list.empty() and !collectionOfLists.empty()) {
		fprintf(stderr, "Options --query/-q and --list/-l and --collectionOfLists/-c cannot be set at the same time.\n");
		exit(EXIT_FAILURE);
	}
	if (!query.empty() and !list.empty()) {
		fprintf(stderr, "Options --query/-q and --list/-l cannot be set at the same time.\n");
		exit(EXIT_FAILURE);
	}
	if (!query.empty() and !collectionOfLists.empty()) {
		fprintf(stderr, "Options --query/-q and --collectionOfLists/-c cannot be set at the same time.\n");
		exit(EXIT_FAILURE);
	}
	if (!list.empty() and !collectionOfLists.empty()) {
		fprintf(stderr, "Options --list/-l and --collectionOfLists/-c cannot be set at the same time.\n");
		exit(EXIT_FAILURE);
	}

	try {
		BarcodesIndex = loadBarcodesIndex(indexFile);
		vector<string> reads;
		if (!query.empty()) {
			if (!gzipped) {
				reads = retrieveReadsWithBarcode(fastqFile, BarcodesIndex, query);
			} else {
				reads = retrieveReadsWithBarcode_Gzip(fastqFile, BarcodesIndex, query);
			}
		} else if (!list.empty()) {
			if (!gzipped) {
				reads = retrieveReadsWithBarcodes(fastqFile, BarcodesIndex, list, nbThreads);
			} else {
				reads = retrieveReadsWithBarcodes_Gzip(fastqFile, BarcodesIndex, list, nbThreads);
			}
		} else {
			in.open(collectionOfLists, ios::in);
			if (!in.is_open()) {
				fprintf(stderr, "Unable to open file %s.", collectionOfLists.c_str());
				exit(EXIT_FAILURE);
			}
			// each line contains a file name containing a list of barcodes (e.g. in this case: listFile = line).
			string listFile;
			while (getline(in, listFile)) {
				if (!gzipped) {
					reads = retrieveReadsWithBarcodes(fastqFile, BarcodesIndex, listFile, nbThreads);
				} else{
					reads = retrieveReadsWithBarcodes_Gzip(fastqFile, BarcodesIndex, listFile, nbThreads);
				}
				// Write the list of reads (variable 'reads') to the corresponding output file (variable 'readsFile').
				string readsFile;
				readsFile = listFile + ".fastq";
				out.open(readsFile, ios::out);
				if (!out.is_open()) {
					fprintf(stderr, "Unable to open file %s.", readsFile.c_str());
					exit(EXIT_FAILURE);
				}
				for (string r : reads) {
					out << r << endl;
				}
				out.close();

				// Est-ce que collectionOfLists correspond en fait a un repertoire avec les listes de barcodes a query dans ce repertoire \
				(et donc pas besoin d'ajouter une etape de creation de fichier)
				// Ou collectionOfLists correspond a un fichier contenant les fichiers de listes de barcodes qu'on veut query \
				(et donc add etape ou on cree ce fichier contenant les fichiers de listes de barcodes)

				// Si dans MTG-Link, on fait le multiprocessing que a partir du gap-filling, comment multiprocesser etapes LRez \
				est-ce que c'est l'etape d'index qui prend du temps, et donc si je charge l'index une seule fois ca suffit \
				ou est-ce que faire la partie "LRez extract" et la partie "LRez query fastq" en sequentiel ca prend du temps
				// Multiprocess LRez avec l'option --threads si necessaire
			}
		}

		// If --query/-q or --list/-l provided, write the list of reads (variable 'reads') to the output file provided by the user or to stdout if no output file is provided.
		if (!query.empty() or !list.empty()) { 
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
	} catch (exception const& e) {
		cerr << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}