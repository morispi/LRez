#include "subcommands/stats.h"
#include "subcommands/help.h"
#include "barcodesExtraction.h"
#include "computeStats.h"
#include <getopt.h>

void subcommandStats(int argc, char* argv[]) {
	string bamFile;
	unsigned numberOfRegions = 1000;
	unsigned regionSize = 1000;
	string outputFile;
	ofstream out;
	unsigned nbThreads = 1;

	const struct option longopts[] = {
		{"bam",				required_argument,	0, 'b'},
		{"regions",			required_argument,	0, 'r'},
		{"size",			required_argument,	0, 's'},
		{"output",			required_argument,	0, 'o'},
		{"threads",			required_argument,	0, 't'},
		{0, 0, 0, 0},
	};
	int index;
	int iarg = 0;

	iarg = getopt_long(argc, argv, "b:r:s:o:t:", longopts, &index);
	if (iarg == -1) {
		subcommandHelp("stats");
	}
	while (iarg != -1) {
		switch (iarg) {
			case 'b':
				bamFile = optarg;
				break;
			case 'r':
				numberOfRegions = stoul(optarg);
				break;
			case 's':
				regionSize = stoul(optarg);
				break;
			case 'o':
				outputFile = optarg;
				break;
			case 't':
				nbThreads = static_cast<uint32_t>(stoul(optarg));
				break;
			default:
				subcommandHelp("stats");
				break;
		}
		iarg = getopt_long(argc, argv, "b:r:s:o:t:", longopts, &index);
	}

	if (bamFile.empty()) {
		fprintf(stderr, "Please specify a BAM file with option --bam/-b.\n");
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

		vector<string> regionsList = extractRegionsList(reader, regionSize);

		reader.Close();

		pair<vector<unsigned>, vector<unsigned>> res = extractBarcodesAndCommonBarcodesCounts(bamFile, regionsList, min(numberOfRegions, (unsigned) regionsList.size()), nbThreads);
		vector<unsigned> barcodesPerRegion = res.first;
		vector<unsigned> commonBarcodes = res.second;

		Stats stats = extractGlobalStats(bamFile, nbThreads);

		if (!outputFile.empty()) {
			out.open(outputFile, ios::out);
			if (!out.is_open()) {
				fprintf(stderr, "Unable to open file %s.", outputFile.c_str());
				exit(EXIT_FAILURE);
			}

			out << "Number of barcodes: " << stats.nbBarcodes << endl;
			out << "Number of mapped reads: " << stats.nbMappedReads << endl;

			out << endl;

			out << "Number of reads per barcode:" << endl;
			out << "\t min: " << stats.readsPerBarcode[0] << endl;
			out << "\t 1st quantile: " << stats.readsPerBarcode[0.25 * stats.readsPerBarcode.size()] << endl;
			out << "\t median: " << stats.readsPerBarcode[0.5 * stats.readsPerBarcode.size()] << endl;
			out << "\t 3rd quantile: " << stats.readsPerBarcode[0.75 * stats.readsPerBarcode.size()] << endl;
			out << "\t max: " << stats.readsPerBarcode[stats.readsPerBarcode.size()-1] << endl;

			out << endl;

			out << "Number of barcodes per region of size " << regionSize << ":" << endl;
			out << "\t min: " << barcodesPerRegion[0] << endl;
			out << "\t 1st quantile: " << barcodesPerRegion[0.25 * barcodesPerRegion.size()] << endl;
			out << "\t median: " << barcodesPerRegion[0.5 * barcodesPerRegion.size()] << endl;
			out << "\t 3rd quantile: " << barcodesPerRegion[0.75 * barcodesPerRegion.size()] << endl;
			out << "\t max: " << barcodesPerRegion[barcodesPerRegion.size()-1] << endl;

			out << endl;

			out << "Number of common barcodes between adjacent regions of size " << regionSize << ":" << endl;
			out << "\t min: " << commonBarcodes[0] << endl;
			out << "\t 1st quantile: " << commonBarcodes[0.25 * commonBarcodes.size()] << endl;
			out << "\t median: " << commonBarcodes[0.5 * commonBarcodes.size()] << endl;
			out << "\t 3rd quantile: " << commonBarcodes[0.75 * commonBarcodes.size()] << endl;
                        out << "\t max: " << commonBarcodes[commonBarcodes.size()-1] << endl;

			out.close();
		} else {
			cout << "Number of barcodes: " << stats.nbBarcodes << endl;
			cout << "Number of mapped reads: " << stats.nbMappedReads << endl;

			cout << endl;

			cout << "Number of reads per barcode:" << endl;
			cout << "\t min: " << stats.readsPerBarcode[0] << endl;
			cout << "\t 1st quantile: " << stats.readsPerBarcode[0.25 * stats.readsPerBarcode.size()] << endl;
			cout << "\t median: " << stats.readsPerBarcode[0.5 * stats.readsPerBarcode.size()] << endl;
			cout << "\t 3rd quantile: " << stats.readsPerBarcode[0.75 * stats.readsPerBarcode.size()] << endl;
			cout << "\t max: " << stats.readsPerBarcode[stats.readsPerBarcode.size()-1] << endl;
			cout << endl;

			cout << "Number of barcodes per region of size " << regionSize << ":" << endl;
                        cout << "\t min: " << barcodesPerRegion[0] << endl;
			cout << "\t 1st quantile: " << barcodesPerRegion[0.25 * barcodesPerRegion.size()] << endl;
			cout << "\t median: " << barcodesPerRegion[0.5 * barcodesPerRegion.size()] << endl;
			cout << "\t 3rd quantile: " << barcodesPerRegion[0.75 * barcodesPerRegion.size()] << endl;
                        cout << "\t max: " << barcodesPerRegion[barcodesPerRegion.size()-1] << endl;

			cout << endl;

			cout << "Number of common barcodes between adjacent regions of size " << regionSize << ":" << endl;
                        cout << "\t min: " << commonBarcodes[0] << endl;
			cout << "\t 1st quantile: " << commonBarcodes[0.25 * commonBarcodes.size()] << endl;
			cout << "\t median: " << commonBarcodes[0.5 * commonBarcodes.size()] << endl;
			cout << "\t 3rd quantile: " << commonBarcodes[0.75 * commonBarcodes.size()] << endl;
                        cout << "\t max: " << commonBarcodes[commonBarcodes.size()-1] << endl;

		}
	} catch (exception const& e) {
		cerr << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}
