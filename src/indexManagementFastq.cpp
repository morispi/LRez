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
		return t[0];
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

	BarcodesIndex barcodesIndex;
	barcode b;
	string barcode, line;

	int64_t pos = in.tellg();
	while (getline(in, line)) {
		barcode = getBXTag(line);
		getline(in, line);
		getline(in, line);
		getline(in, line);
		if (isValidBarcode(barcode)) {
			b = stringToBarcode(barcode);
			barcodesIndex[b].push_back(pos);
		}
		pos = in.tellg();
	}

	in.close();
	return barcodesIndex;
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

    BarcodesIndex barcodesIndex;
	barcode b;
	string barcode, line;

	int64_t pos = 0;
	int64_t oldPos = 0;

	while (getline(in, line)) {
		barcode = getBXTag(line);
		pos += line.length() + 1;
		getline(in, line);
		pos += line.length()  + 1;
		getline(in, line);
		pos += line.length()  + 1;
		getline(in, line);
		pos += line.length()  + 1;
		if (isValidBarcode(barcode)) {
			b = stringToBarcode(barcode);
			barcodesIndex[b].push_back(oldPos);
		}
		oldPos = pos;
	}

	file.close();
	return barcodesIndex;
}


void saveBarcodesIndex(BarcodesIndex& barcodesIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		fprintf(stderr, "Unable to open file %s.", file.c_str());
		exit(EXIT_FAILURE);
	}

	// Write number of bits per barcode
	unsigned bitsPerBarcode = barcodesIndex.begin()->first.size();
	out << bitsPerBarcode << endl;

	for (auto i : barcodesIndex) {
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