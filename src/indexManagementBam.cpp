#include "indexManagementBam.h"
#include "reverseComplement.h"
#include "../CTPL/ctpl_stl.h"

/**
	Index the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
BarcodesOffsetsIndex indexBarcodesOffsetsFromBamRegions(int id, BamReader& reader, bool primary, unsigned quality, string begRegion, string endRegion, bool tilEnd = 0) {
	BarcodesOffsetsIndex barcodesOffsetsIndex;
	BamAlignment al;
	barcode b;
	string bc, barcode;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region" + region + ".");
	}

	// Process and store alignments located between the two regions of interest
	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment and (!primary or (primary and al.IsMapped() and al.IsPrimaryAlignment())) and (quality == 0 or (al.IsMapped() and al.MapQuality >= quality))) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodesOffsetsIndex[b].push_back(pos);
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
			if ((!primary or (primary and al.IsMapped() and al.IsPrimaryAlignment())) and (quality == 0 or (al.IsMapped() and al.MapQuality >= quality))) {
				barcode.clear();
				al.GetTag(BXTAG, barcode);
				if (isValidBarcode(barcode)) {
					b = stringToBarcode(barcode);
					barcodesOffsetsIndex[b].push_back(pos);
				}
			}
			pos = reader.FTell();
		}
	}

	return barcodesOffsetsIndex;
}

BarcodesOffsetsIndex indexBarcodesOffsetsFromBam(string bamFile, bool primary, unsigned quality, unsigned nbThreads) {
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
	vector<std::future<BarcodesOffsetsIndex>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(indexBarcodesOffsetsFromBamRegions, ref(readers[i]), ref(primary), ref(quality), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	BarcodesOffsetsIndex barcodesOffsetsIndex;
	BarcodesOffsetsIndex curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			for (auto v : p.second) {
				barcodesOffsetsIndex[p.first].push_back(v);
			}
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return barcodesOffsetsIndex;
}


void saveBarcodesOffsetsIndex(BarcodesOffsetsIndex& barcodesOffsetsIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		throw ios_base::failure("open: Unable to open file " + file + ". Please make sure the files exists.");
	}

	// Write number of bits per barcode
	unsigned bitsPerBarcode = barcodesOffsetsIndex.begin()->first.size();
	out << bitsPerBarcode << endl;

	for (auto i : barcodesOffsetsIndex) {
		// Convert bool vector to int value, and write barcode
		uint64_t barcode = 0ULL;
		uint64_t tmp;
		for (unsigned j = 0; j < bitsPerBarcode; j++) {
			tmp = i.first[j];
			barcode |= tmp << (bitsPerBarcode - j - 1);
		}
		out << barcode << ";";

		// Write BamRegion vector
		for (unsigned j = 0; j < i.second.size() - 1; j++) {
			out << i.second[j] << ","; 
		}
		out << i.second[i.second.size() - 1] << endl;
	}

	out.close();
} 

BarcodesOffsetsIndex loadBarcodesOffsetsIndex(string file) {
	ifstream in;
	in.open(file);
	if (!in.is_open()) {
		throw ios_base::failure("open: Unable to barcodes index file " + file + ". Please provide an existing and valid file.");
	}

	BarcodesOffsetsIndex index;
	string line;
	vector<string> v, w, r;
	barcode barcode;

	// Get number of bits per barcode
	getline(in, line);
	unsigned bitsPerBarcode = stoul(line);

	while (getline(in, line)) {
		// Get barcode
		v = splitString(line, ";");
		uint64_t b = stoull(v[0]);
		barcode.clear();
		for (unsigned i = 0; i < bitsPerBarcode; i++) {
			barcode.push_back(b & (1ULL << (bitsPerBarcode - i - 1)));
		}

		// Get offsets
		w = splitString(v[1], ",");
		vector<int64_t> regions;
		for (string ww : w) {
			regions.push_back(stoul(ww));
		}
		index[barcode] = regions;
	}

	in.close();
	return index;
}

/**
	Index the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
BarcodesPositionsIndex indexBarcodesPositionsFromBamRegions(int id, BamReader& reader, bool primary, unsigned quality, string begRegion, string endRegion, bool tilEnd = 0) {
	BarcodesPositionsIndex barcodesPositionsIndex;
	BamAlignment al;
	barcode b;
	string bc, barcode;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region" + region + ".");
	}

	// Process and store alignments located between the two regions of interest
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = al.IsMapped() and (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment and (!primary or (primary and al.IsMapped() and al.IsPrimaryAlignment())) and (quality == 0 or (al.IsMapped() and al.MapQuality >= quality))) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodesPositionsIndex[b].push_back(make_pair(al.RefID, al.Position));
			}
		}
	}

	return barcodesPositionsIndex;
}


BarcodesPositionsIndex indexBarcodesPositionsFromBam(string bamFile, bool primary, unsigned quality, unsigned nbThreads) {
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

	// Split the processing of the BAM files into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<BarcodesPositionsIndex>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(indexBarcodesPositionsFromBamRegions, ref(readers[i]), ref(primary), ref(quality), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	BarcodesPositionsIndex barcodesPositionsIndex;
	BarcodesPositionsIndex curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			for (auto v : p.second) {
				barcodesPositionsIndex[p.first].push_back(v);
			}
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return barcodesPositionsIndex;
}

void saveBarcodesPositionsIndex(BarcodesPositionsIndex& barcodesPositionsIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		throw ios_base::failure("open: Unable to open file " + file + ". Please make sure the files exists.");
	}

	// Write number of bits per barcode
	unsigned bitsPerBarcode = barcodesPositionsIndex.begin()->first.size();
	out << bitsPerBarcode << endl;

	for (auto i : barcodesPositionsIndex) {
		// Convert bool vector to int value, and write barcode
		uint64_t barcode = 0ULL;
		uint64_t tmp;
		for (unsigned j = 0; j < bitsPerBarcode; j++) {
			tmp = i.first[j];
			barcode |= tmp << (bitsPerBarcode - j - 1);
		}
		out << barcode << ";";

		// Write BamRegion vector
		for (unsigned j = 0; j < i.second.size() - 1; j++) {
			out << i.second[j].first << ":" << i.second[j].second << ","; 
		}
		out << i.second[i.second.size() - 1].first << ":" << i.second[i.second.size() - 1].second << endl;
	}

	out.close();
}

BarcodesPositionsIndex loadBarcodesPositionsIndex(string file) {
	ifstream in;
	in.open(file);
	if (!in.is_open()) {
		throw ios_base::failure("open: Unable to barcodes index file " + file + ". Please provide an existing and valid file.");
	}

	BarcodesPositionsIndex index;
	string line;
	vector<string> v, w, r;
	barcode barcode;

	// Get number of bits per barcode
	getline(in, line);
	unsigned bitsPerBarcode = stoul(line);

	while (getline(in, line)) {
		// Get barcode
		v = splitString(line, ";");
		uint64_t b = stoull(v[0]);
		barcode.clear();
		for (unsigned i = 0; i < bitsPerBarcode; i++) {
			barcode.push_back(b & (1ULL << (bitsPerBarcode - i - 1)));
		}

		// Get positions
		w = splitString(v[1], ",");
		vector<pair<int32_t, int32_t>> positions;
		for (string ww : w) {
			vector<string> p = splitString(ww, ":");
			positions.push_back(make_pair(stoi(p[0]), stoi(p[1])));
		}
		index[barcode] = positions;
	}

	in.close();
	return index;
}
