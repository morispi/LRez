#include "barcodesComparison.h"
#include "barcodesExtraction.h"
#include "alignmentsRetrieval.h"

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
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	ifstream in;
	in.open(regions);
	if (!in.is_open()) {
		fprintf(stderr, "Unable to open barcodes index file %s. Please provide an existing and valid file.\n", regions.c_str());
		exit(EXIT_FAILURE);
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
			// TODO : peut être une erreur ici ? ne comparerait qu'une des extrémities d'un contig cible si un même barcode est présent aux deux extrémités
			// Plutôt faire un set de visited avec les tRegion ?
			if (rv[al.RefID].RefLength < size) {
				tRegion = rv[al.RefID].RefName + ":0-" + rv[al.RefID].RefName + ":" + to_string(rv[al.RefID].RefLength);
				if (!consideredRegions.count(tRegion)) {
					counts[make_pair(qRegion, tRegion)]++;
					consideredRegions.insert(tRegion);
				}
			} else if (al.Position < size) {
				tRegion = rv[al.RefID].RefName + ":0-" + rv[al.RefID].RefName + ":" + to_string(size);
				if (!consideredRegions.count(tRegion)) {
					counts[make_pair(qRegion, tRegion)]++;
					consideredRegions.insert(tRegion);
				}
			} else if (rv[al.RefID].RefLength - size < al.Position) {
				tRegion = rv[al.RefID].RefName + ":" + to_string(rv[al.RefID].RefLength - size) + "-" + rv[al.RefID].RefName + ":" + to_string(rv[al.RefID].RefLength);
				if (!consideredRegions.count(tRegion)) {
					counts[make_pair(qRegion, tRegion)]++;
					consideredRegions.insert(tRegion);
				}
			}
		}
	}
}

robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContig(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contig, int size) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	int id = reader.GetReferenceID(contig);
	if (id == -1) {
		fprintf(stderr, "Cannot find refence with ID %s.\n", contig.c_str());
		exit(EXIT_FAILURE);
	}

	robin_hood::unordered_set<barcode> barcodes;
	robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> counts;
	const RefVector rv = reader.GetReferenceData();
	string qRegion;
	
	// Process left extremity of the contig
	if (rv[id].RefLength < size) {
		qRegion = contig + ":0-" + contig + ":" + to_string(rv[id].RefLength);
	} else {
		qRegion = contig + ":0-" + contig + ":" + to_string(size);
	}
	computeCommonBarcodesCounts(counts, BarcodesOffsetsIndex, reader, rv, size, qRegion);

	// Process left extremity of the contig, if it is long enough
	if (rv[id].RefLength > size) {
		qRegion = contig + ":" + to_string(rv[id].RefLength - size) + "-" + contig + ":" + to_string(rv[id].RefLength);
		computeCommonBarcodesCounts(counts, BarcodesOffsetsIndex, reader, rv, size, qRegion);
	}

	reader.Close();
	return counts;

}