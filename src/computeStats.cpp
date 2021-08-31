#include "computeStats.h"
#include "barcodesExtraction.h"
#include "barcodesComparison.h"
#include "../CTPL/ctpl_stl.h"

vector<unsigned> extractBarcodesCountsPerRegion(BamReader& reader, vector<string> regionsList, unsigned numberOfRegions) {
	vector<unsigned> res(numberOfRegions);

	srand(0);
	for (unsigned i = 0; i < numberOfRegions; i++) {
		unsigned n = rand() % (numberOfRegions - 2) + 1;
		res[i] = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n]).size();
	}

	std::sort(res.begin(), res.end());
	return res;
}

vector<unsigned> extractCommonBarcodesCounts(BamReader& reader, vector<string> regionsList, unsigned numberOfRegions) {
	vector<unsigned> res(2 * numberOfRegions);

	srand(0);
	for (unsigned i = 0; i < numberOfRegions; i++) {
		unsigned n = rand() % (numberOfRegions - 2) + 1;
		robin_hood::unordered_set<barcode> barcodesRegion = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n]);
		robin_hood::unordered_set<barcode> barcodesPreviousRegion = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n-1]);
		robin_hood::unordered_set<barcode> barcodesFollowingRegion = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n+1]);
		res[2 * i] = countCommonBarcodes(barcodesRegion, barcodesPreviousRegion);
		res[2 * i + 1] = countCommonBarcodes(barcodesRegion, barcodesFollowingRegion);
	}

	std::sort(res.begin(), res.end());
	return res;
}

pair<vector<unsigned>, vector<unsigned>> extractBarcodesAndCommonBarcodesCountsFromRegionsList(int id, BamReader& reader, vector<string>& regionsList, vector<unsigned>& regionsIds) {
	vector<unsigned> counts;
	vector<unsigned> common;

	for (unsigned n : regionsIds) {
		robin_hood::unordered_set<barcode> barcodesRegion = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n]);
		robin_hood::unordered_set<barcode> barcodesPreviousRegion = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n-1]);
		robin_hood::unordered_set<barcode> barcodesFollowingRegion = extractBarcodesBitsFromRegion_BamReader(reader, regionsList[n+1]);

		counts.push_back(barcodesRegion.size());
		common.push_back(countCommonBarcodes(barcodesRegion, barcodesPreviousRegion));
		common.push_back(countCommonBarcodes(barcodesRegion, barcodesFollowingRegion));
	}

	return make_pair(counts, common);
}

pair<vector<unsigned>, vector<unsigned>> extractBarcodesAndCommonBarcodesCounts(string bamFile, vector<string> regionsList, unsigned numberOfRegions, unsigned nbThreads) {
	// Open the necessary number of BamReaders
	vector<BamReader> readers(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		if (!readers[i].Open(bamFile)) {
			throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
		}
		if (!readers[i].LocateIndex()) {
			throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
		}
		readers[i].Rewind();
	}

	// Fill vectors with the ids of the regions to process
	vector<vector<unsigned>> queries(nbThreads);
	unsigned curThread = 0;
	string line;
	for (unsigned i = 0; i < numberOfRegions; i++) {
		unsigned n = rand() % (numberOfRegions - 2) + 1;
		queries[curThread].push_back(n);
		curThread = (curThread + 1) % nbThreads;
    }

    // Split the processing into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<pair<vector<unsigned>, vector<unsigned>>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(extractBarcodesAndCommonBarcodesCountsFromRegionsList, ref(readers[i]), ref(regionsList), ref(queries[i]));
	}

	// Retrieve threads subresults and build global result
	vector<unsigned> counts;
	vector<unsigned> common;
	pair<vector<unsigned>, vector<unsigned>> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		
		for (auto c : curRes.first) {
			counts.push_back(c);
		}

		for (auto c : curRes.second) {
			common.push_back(c);
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	std::sort(counts.begin(), counts.end());
	std::sort(common.begin(), common.end());
	return make_pair(counts, common);
}

// Structure to store temporary stats
struct AuxStats {
	unsigned nbMappedReads;
	robin_hood::unordered_map<barcode, robin_hood::unordered_set<string>> barcodeToReadSet;
};

/**
	Compute stats from the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
AuxStats extractGlobalStatsFromBamRegions(int id, BamReader& reader, string begRegion, string endRegion, bool tilEnd = 0) {
	AuxStats stats;
	stats.nbMappedReads = 0;
	BamAlignment al;
	barcode b;
	string bc;
	robin_hood::unordered_set<string> mappedReads;
	robin_hood::unordered_map<barcode, robin_hood::unordered_set<string>> barcodeToReadSet;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	reader.SetRegion(r);

	// Process and extract barcodes from alignments located between the two regions of interest
	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment) {
			bc.clear();
			al.GetTag(BXTAG, bc);
			if (isValidBarcode(bc)) {
				b = stringToBarcode(bc);
				stats.barcodeToReadSet[b].insert(al.Name);
			}

			if (al.IsMapped() and al.IsPrimaryAlignment() and !(al.AlignmentFlag & Constants::BAM_ALIGNMENT_SUPPLEMENTARY)) {
				stats.nbMappedReads += 1;
			}
		}
		pos = reader.FTell();
	}

	// Treat unmapped reads if we're at the end of the file
	if (tilEnd) {
		reader.Rewind();
		if (!reader.FSeek(pos)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(pos) + ".");
		}
		while (reader.GetNextAlignment(al)) {
			bc.clear();
			al.GetTag(BXTAG, bc);
			if (isValidBarcode(bc)) {
				b = stringToBarcode(bc);
				stats.barcodeToReadSet[b].insert(al.Name);
			}

			if (al.IsMapped() and al.IsPrimaryAlignment() and !(al.AlignmentFlag & Constants::BAM_ALIGNMENT_SUPPLEMENTARY)) {
				stats.nbMappedReads += 1;
			}
			pos = reader.FTell();
		}
	}

	return stats;
}

Stats extractGlobalStats(string bamFile, unsigned nbThreads) {
	// Open the necessary number of BamReaders
	vector<BamReader> readers(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		if (!readers[i].Open(bamFile)) {
			throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
		}
		if (!readers[i].LocateIndex()) {
			throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
		}
		readers[i].Rewind();
	}

	// Split the reference genome into regions of size 10000
	vector<string> regions = extractRegionsList(readers[0], 10000);
	readers[0].Rewind();

	// Split the processing of the BAM file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<AuxStats>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(extractGlobalStatsFromBamRegions, ref(readers[i]), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	Stats stats;
	stats.nbBarcodes = 0;
	robin_hood::unordered_map<barcode, robin_hood::unordered_set<string>> barcodeToReadSet;
	AuxStats curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		stats.nbMappedReads += curRes.nbMappedReads;

		for (auto p : curRes.barcodeToReadSet) {
			for (auto pp : p.second) {
				barcodeToReadSet[p.first].insert(pp);
			}
		}
	}

	stats.nbBarcodes = barcodeToReadSet.size();

	for (auto p : barcodeToReadSet) {
		stats.readsPerBarcode.push_back(p.second.size());
	}
	std::sort(stats.readsPerBarcode.begin(), stats.readsPerBarcode.end());


	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return stats;
}