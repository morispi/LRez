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

	BarcodesOffsetsIndex BarcodesOffsetsIndex;
	BamAlignment al;
	barcode b;
	string bc, barcode;
	BamRegion r;

	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		if (al.IsMapped() and (!primary or (primary and al.IsPrimaryAlignment())) and al.MapQuality >= quality) {
			barcode.clear();
			// LongRanger
			if (al.GetTag(BXTAG, bc)) {
				barcode = bc.substr(0, BARCODE_SIZE);
			} else if (CONSIDER_RX and al.GetTag(RXTAG, bc)) {
				barcode = bc.substr(0, BARCODE_SIZE);
			// Other aligners
			} else if (!al.GetTag(BXTAG, bc) and !al.GetTag(RXTAG, bc)) {
				// Barcode is only contained in first mate, so only index if first mate or unpaired
				if (!al.IsPaired() or al.IsFirstMate()) {
					if (!al.IsReverseStrand()) {
						barcode = al.QueryBases.substr(0, BARCODE_SIZE);
					} else {
						barcode = rev_comp::run(al.QueryBases).substr(0, BARCODE_SIZE);
					}				
				}
			}
			if (!barcode.empty()) {
				b = stringToBarcode(barcode);
				BarcodesOffsetsIndex[b].push_back(pos);
			}
		}
		pos = reader.FTell();
	}

	reader.Close();
	return BarcodesOffsetsIndex;
}


void saveBarcodesOffsetsIndex(BarcodesOffsetsIndex& BarcodesOffsetsIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", file.c_str());
		exit(EXIT_FAILURE);
	}

	for (auto i : BarcodesOffsetsIndex) {
		// Write barcode
		out << i.first << ";";

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

	while (getline(in, line)) {
		v = splitString(line, ";");
		barcode = stol(v[0]);
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
			// LongRanger
			if (al.GetTag(BXTAG, bc)) {
				barcode = bc.substr(0, BARCODE_SIZE);
			} else if (CONSIDER_RX and al.GetTag(RXTAG, bc)) {
				barcode = bc.substr(0, BARCODE_SIZE);
			// Other aligners
			} else if (!al.GetTag(BXTAG, bc) and !al.GetTag(RXTAG, bc)) {
				// Barcode is only contained in first mate, so only index if first mate or unpaired
				if (!al.IsPaired() or al.IsFirstMate()) {
					if (!al.IsReverseStrand()) {
						barcode = al.QueryBases.substr(0, BARCODE_SIZE);
					} else {
						barcode = rev_comp::run(al.QueryBases).substr(0, BARCODE_SIZE);
					}				
				}
			}
			if (!barcode.empty()) {
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

	for (auto i : barcodesPositionsIndex) {
		// Write barcode
		// out << barcodeToString(i.first, BARCODE_SIZE) << ";";
		out << i.first << ";";

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

	while (getline(in, line)) {
		v = splitString(line, ";");
		// barcode = stringToBarcode(v[0]);
		barcode = stol(v[0]);
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
