#include "indexManagementBam.h"
#include "reverseComplement.h"

BarcodesOffsetsIndex indexBarcodesOffsetsFromBam(string bamFile, bool primary, unsigned quality) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	BarcodesOffsetsIndex barcodesOffsetsIndex;
	BamAlignment al;
	barcode b;
	string bc, barcode;
	BamRegion r;

	int64_t pos = reader.FTell();
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

	reader.Close();
	return barcodesOffsetsIndex;
}


void saveBarcodesOffsetsIndex(BarcodesOffsetsIndex& barcodesOffsetsIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", file.c_str());
		exit(EXIT_FAILURE);
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
		fprintf(stderr, "Unable to open barcodes index file %s. Please provide an existing and valid file.\n", file.c_str());
		exit(EXIT_FAILURE);
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

BarcodesPositionsIndex indexBarcodesPositionsFromBam(string bamFile, bool primary, unsigned quality) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	BarcodesPositionsIndex barcodesPositionsIndex;
	BamAlignment al;
	barcode b;
	string bc, barcode;
	BamRegion r;

	while (reader.GetNextAlignment(al)) {
		if (al.IsMapped() and (!primary or (primary and al.IsPrimaryAlignment())) and al.MapQuality >= quality) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodesPositionsIndex[b].push_back(make_pair(al.RefID, al.Position));
			}
		}
	}

	reader.Close();
	return barcodesPositionsIndex;
}

void saveBarcodesPositionsIndex(BarcodesPositionsIndex& barcodesPositionsIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", file.c_str());
		exit(EXIT_FAILURE);
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
		fprintf(stderr, "Unable to open barcodes index file %s. Please provide an existing and valid file.\n", file.c_str());
		exit(EXIT_FAILURE);
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
