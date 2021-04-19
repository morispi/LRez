#include "indexManagementFastq.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

string getBXTag(string header) {
	vector<string> t = splitString(header, "BX:Z:");
	t = splitString(t[1], "\t");
	t = splitString(t[0], " ");

	if (t.size() == 1 and t[0].length() != 0) {
		return splitString(t[0], "-1")[0];
	} else {
		return "";
	}
}

string getRXTag(string header) {
	vector<string> t = splitString(header, "RX:Z:");
	t = splitString(t[1], "\t");
	t = splitString(t[0], " ");

	if (t.size() == 1 and t[0].length() != 0) {
		return splitString(t[0], "-1")[0];
	} else {
		return "";
	}
}

BarcodesIndex indexBarcodesFromFastq(string fastqFile) {
	ifstream in;
	in.open(fastqFile);
	if (!in.is_open()) {
		fprintf(stderr, "Unable to open fastq file %s. Please provide an existing and valid file.\n", fastqFile.c_str());
		exit(EXIT_FAILURE);
	}

	BarcodesIndex BarcodesIndex;
	barcode b;
	string barcode, line;

	int64_t pos = in.tellg();
	while (getline(in, line)) {
		barcode = getBXTag(line);
		if (barcode.empty() and CONSIDER_RX) {
			barcode = getRXTag(line);
		}
		getline(in, line);
		getline(in, line);
		getline(in, line);
		if (!barcode.empty()) {
			b = stringToBarcode(barcode);
			BarcodesIndex[b].push_back(pos);
		}
		pos = in.tellg();
	}

	in.close();
	return BarcodesIndex;
}

BarcodesIndex indexBarcodesFromFastqGz(string fastqFile) {
	std::ifstream file(fastqFile, std::ios_base::in | std::ios_base::binary);
	if (!file.is_open()) {
		fprintf(stderr, "Unable to open fastq file %s. Please provide an existing and valid file.\n", fastqFile.c_str());
		exit(EXIT_FAILURE);
	}

    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    istream in(&inbuf);

    BarcodesIndex BarcodesIndex;
	barcode b;
	string barcode, line;

	// int64_t pos = in.tellg();
	int64_t pos = 0;
	int64_t oldPos = 0;

	while (getline(in, line)) {
		barcode = getBXTag(line);
		pos += line.length() + 1;
		if (barcode.empty() and CONSIDER_RX) {
			barcode = getRXTag(line);
		}
		getline(in, line);
		pos += line.length()  + 1;
		getline(in, line);
		pos += line.length()  + 1;
		getline(in, line);
		pos += line.length()  + 1;
		if (!barcode.empty()) {
			b = stringToBarcode(barcode);
			BarcodesIndex[b].push_back(oldPos);
		}
		oldPos = pos;
	}

	file.close();
	return BarcodesIndex;
}


void saveBarcodesIndex(BarcodesIndex& BarcodesIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", file.c_str());
		exit(EXIT_FAILURE);
	}

	for (auto i : BarcodesIndex) {
		// Write barcode
		// out << barcodeToString(i.first, BARCODE_SIZE) << ";";
		// out << i.first << ";";
		// for (bool b : i.first) {
		// 	out << b;
		// }
		// out << ";";
		out << i.first << ";";

		// Write BamRegion vector
		for (unsigned j = 0; j < i.second.size() - 1; j++) {
			out << i.second[j] << ","; 
		}
		out << i.second[i.second.size() - 1] << endl;
	}

	out.close();
}

BarcodesIndex loadBarcodesIndex(string file) {
	ifstream in;
	in.open(file);
	if (!in.is_open()) {
		fprintf(stderr, "Unable to open barcodes index file %s. Please provide an existing and valid file.\n", file.c_str());
		exit(EXIT_FAILURE);
	}

	BarcodesIndex index;
	string line;
	vector<string> v, w, r;
	barcode barcode;

	while (getline(in, line)) {
		v = splitString(line, ";");
		// barcode = stringToBarcode(v[0]);
		// barcode = stol(v[0]);
		// barcode.clear();
		// for (char c : v[0]) {
		// 	if (c == '0') {
		// 		barcode.push_back(false);
		// 	} else {
		// 		barcode.push_back(true);
		// 	}
		// }
		barcode = v[0];
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