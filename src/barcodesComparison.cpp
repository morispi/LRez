#include "barcodesComparison.h"
#include "barcodesExtraction.h"
#include "alignmentsRetrieval.h"
#include "../CTPL/ctpl_stl.h"

unsigned countCommonBarcodes(robin_hood::unordered_set<barcode> barcodes1, robin_hood::unordered_set<barcode> barcodes2) {
	unsigned res = 0;

	for (const auto b : barcodes1) {
		 res += barcodes2.count(b);
    }

	return res;
}

robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> computePairwiseCommonBarcounts(robin_hood::unordered_map<string, robin_hood::unordered_set<barcode>> regionsBarcodes) {
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> res;
	
	auto it1 = regionsBarcodes.begin();
	while (it1 != regionsBarcodes.end()) {
		auto it2 = it1;
		it2++;
		while (it2 != regionsBarcodes.end()) {
			res[make_pair(it1->first, it2->first)] = countCommonBarcodes(it1->second, it2->second);
			it2++;
		}
		it1++;
	}

	return res;
}

robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareRegions(string bamFile, string regions) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
	}

	ifstream in;
	in.open(regions);
	if (!in.is_open()) {
		throw ios_base::failure("open: Unable to open regions file " + regions + ". Please provide an existing and valid file.");
	}

	// Fill map with barcodes of each region
	robin_hood::unordered_map<string, robin_hood::unordered_set<barcode>> regionsBarcodes;
	string region;
	while (getline(in, region)) {
		regionsBarcodes[region] = extractBarcodesBitsFromRegion_BamReader(reader, region);
	}

	reader.Close();
	in.close();

	// Count common barcodes of each pair of regions
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> res = computePairwiseCommonBarcounts(regionsBarcodes);

	return res;
}

void computeCommonBarcodesCounts(robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs>& counts, BarcodesOffsetsIndex& BarcodesOffsetsIndex, BamReader& reader, const RefVector& rv, int size, string& qRegion) {
	robin_hood::unordered_set<string> consideredRegions;
	string tRegion;
	robin_hood::unordered_set<barcode> barcodes;

	barcodes = extractBarcodesBitsFromRegion_BamReader(reader, qRegion);
	for (auto b : barcodes) {
		vector<BamAlignment> alignments = retrieveAlignmentsWithBarcodeBits_BamReader(reader, BarcodesOffsetsIndex, b);
		consideredRegions.clear();
		for (auto al : alignments) {
			if (al.IsMapped() and rv[al.RefID].RefLength < size) {
				tRegion = rv[al.RefID].RefName + ":0-" + rv[al.RefID].RefName + ":" + to_string(rv[al.RefID].RefLength);
				if (!consideredRegions.count(tRegion)) {
					counts[make_pair(qRegion, tRegion)]++;
					consideredRegions.insert(tRegion);
				}
			} else if (al.IsMapped() and al.Position < size) {
				tRegion = rv[al.RefID].RefName + ":0-" + rv[al.RefID].RefName + ":" + to_string(size);
				if (!consideredRegions.count(tRegion)) {
					counts[make_pair(qRegion, tRegion)]++;
					consideredRegions.insert(tRegion);
				}
			} else if (al.IsMapped() and rv[al.RefID].RefLength - size < al.Position) {
				tRegion = rv[al.RefID].RefName + ":" + to_string(rv[al.RefID].RefLength - size) + "-" + rv[al.RefID].RefName + ":" + to_string(rv[al.RefID].RefLength);
				if (!consideredRegions.count(tRegion)) {
					counts[make_pair(qRegion, tRegion)]++;
					consideredRegions.insert(tRegion);
				}
			}
		}
	}
}

robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContig_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contig, int size) {
	int id = reader.GetReferenceID(contig);
	if (id == -1) {
		throw runtime_error("GetReferenceID: Cannot find reference with ID " + contig + ".");
	}

	robin_hood::unordered_set<barcode> barcodes;
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> counts;
	const RefVector rv = reader.GetReferenceData();
	string qRegion;
	
	// Process left extremity of the contig
	if (rv[id].RefLength < size) {
		qRegion = contig + ":0-" + to_string(rv[id].RefLength);
	} else {
		qRegion = contig + ":0-" + to_string(size);
	}
	computeCommonBarcodesCounts(counts, BarcodesOffsetsIndex, reader, rv, size, qRegion);

	// Process left extremity of the contig, if it is long enough
	if (rv[id].RefLength > size) {
		qRegion = contig + ":" + to_string(rv[id].RefLength - size) + "-" + to_string(rv[id].RefLength);
		computeCommonBarcodesCounts(counts, BarcodesOffsetsIndex, reader, rv, size, qRegion);
	}

	return counts;
}

robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContig(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contig, int size) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");

	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
	}

	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> res = compareContig_BamReader(reader, BarcodesOffsetsIndex, contig, size);
	reader.Close();
	
	return res;
}

/**
	Compare the extremities of a given list of contigs to all other contigs' extremities
*/
robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContigsList(int id, BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, int size, vector<string>& queries) {
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> res;
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> tmpRes;
	for (string q : queries) {
		tmpRes = compareContig_BamReader(reader, BarcodesOffsetsIndex, q, size);
		for (auto r : tmpRes) {
			res[r.first] += r.second;
		}
	}

	return res;
}

robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContigs(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contigs, int size, unsigned nbThreads) {
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

	// Fill vectors with the query barcodes
	ifstream bc;
	bc.open(contigs);
	if (!bc.is_open()) {
		throw ios_base::failure("open: Unable to open contigs list file " + contigs + ". Please provide an existing and valid file.");
	}

	vector<vector<string>> queries(nbThreads);
	unsigned curThread = 0;
	string line;
	while (getline(bc, line)) {
		queries[curThread].push_back(line);
		curThread = (curThread + 1) % nbThreads;
    }

    bc.close();

	// Split the querying into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(compareContigsList, ref(readers[i]), ref(BarcodesOffsetsIndex), size, ref(queries[i]));
	}

	// Retrieve threads subresults and build global result
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> res;
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto r : curRes) {
			res[r.first] += r.second;
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return res;
}