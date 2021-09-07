#include "indexManagementFastq.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include "../CTPL/ctpl_stl.h"
#include "gzIndex.h"

string getBXTag(string header) {
	vector<string> t = splitString(header, "BX:Z:");
	if (t.size() != 2) {
		return "";
	}

	t = splitString(t[1], " ");
	t = splitString(t[0], "\n");
	t = splitString(t[0], " ");
	t = splitString(t[0], "\n");

	if (t.size() == 1 and t[0].length() != 0) {
		return t[0];
	} else {
		return "";
	}
}

/**
	Index the portion of the FASTQ file located between begOffset and endOffset.
*/
BarcodesIndex indexBarcodesFromFastqOffsets(int id, ifstream& in, streampos begOffset, streampos endOffset) {
	BarcodesIndex barcodesIndex;
	barcode b;
	string barcode, line;

	in.seekg(begOffset, ios_base::beg);
	int64_t pos = in.tellg();

	// Read until the first read to process
	getline(in, line);
	while (getBXTag(line) == "" or line[0] != '@') {
		pos = in.tellg(); 
		getline(in, line);
	}

	// Process the first read
	barcode = getBXTag(line);
	getline(in, line);
	getline(in, line);
	getline(in, line);
	if (isValidBarcode(barcode)) {
		b = stringToBarcode(barcode);
		barcodesIndex[b].push_back(pos);
	}
	pos = in.tellg();

	// Process remaining reads until the end offset
	while (pos <= endOffset and getline(in, line)) {
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

	return barcodesIndex;
}

BarcodesIndex indexBarcodesFromFastq(string fastqFile, unsigned nbThreads) {
	// Open the necessary number of ifstreams
	vector<ifstream> ins(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		ins[i].open(fastqFile);
		if (!ins[i].is_open()) {
			throw ios_base::failure("open: Unable to open fastq file " + fastqFile + ". Please provide an existing and valid file.");
		}
	}

	// Get the max offset of the FASTQ file
	ins[0].seekg(0, ios_base::end);
	streampos maxOffset = ins[0].tellg();
	ins[0].seekg(0, ios_base::beg);

	// Split the processing of the fastq file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<BarcodesIndex>> results(nbThreads);
	int64_t curLimit = 0;
	int64_t newLimit = 0;
	for (unsigned i = 0; i < nbThreads; i++) {
		newLimit = (int64_t) ((float) maxOffset / nbThreads) * (i + 1);
		results[i] = myPool.push(indexBarcodesFromFastqOffsets, ref(ins[i]), curLimit, newLimit);
		curLimit = newLimit + 1;
	}

	// Retrieve threads subresults and build global index
	BarcodesIndex barcodesIndex;
	BarcodesIndex curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			for (auto v : p.second) {
				barcodesIndex[p.first].push_back(v);
			}
		}
	}

	// Close ifstreams
	for (unsigned i = 0; i < nbThreads; i++) {
		ins[i].close();
	}

	return barcodesIndex;
}

/**
	Read the buffer of size size starting from pos, and up to the next newline symbol.
	Reads a new chunk of the gzipped file in and stores it in gzIndex if needed.
	Return the line that was read.
*/
string nextLine(FILE* in, struct access* gzIndex, unsigned char* buffer, unsigned size, off_t& pos, off_t& totalPos) {
	string line = "";

	while (pos < size and buffer[pos] != '\n') {
		line += buffer[pos];
		pos++;
		totalPos++;
	}

	// Read a new chunk of the file if no newline symbol was found
	if (pos == size and line[line.length() - 1] != '\n') {
		extract_Stream(in, gzIndex, totalPos, buffer, BUFSIZE);
		pos = 0;
		
		while (pos < size and buffer[pos] != '\n') {
			line += buffer[pos];
			pos++;
			totalPos++;
		}
	}


	pos++;
	totalPos++;
	return line;
}


/**
	Index the portion of the FASTQ file located between begOffset and endOffset.
	tilEnd has to be set to one for the thread processing the end of the FASTQ file, so it can process it beyond the last window of the index.
*/
BarcodesIndex indexBarcodesFromFastqGzOffsets(int id, FILE* in, struct access* gzIndex, int firstWindow, int lastWindow, bool tilEnd = 0) {
	BarcodesIndex barcodesIndex;
	barcode b;
	string barcode, line;
	off_t curPos = 0;
	off_t oldPos = 0;
	off_t pos = gzIndex->list[firstWindow].out;

	unsigned char* buf = (unsigned char*) malloc(BUFSIZE);
	extract_Stream(in, gzIndex, pos, buf, BUFSIZE);

	// Read until the first read to process
	line = nextLine(in, gzIndex, buf, BUFSIZE, curPos, pos);
	while (getBXTag(line) == "" or line[0] != '@') {
		oldPos = pos;
		line = nextLine(in, gzIndex, buf, BUFSIZE, curPos, pos);
	}

	while ((pos <= gzIndex->list[lastWindow].out or (tilEnd and pos < gzIndex->maxOffset)) and line != "") {
		barcode = getBXTag(line);
		if (isValidBarcode(barcode)) {
			b = stringToBarcode(barcode);
			barcodesIndex[b].push_back(oldPos);
		}
		
		line = nextLine(in, gzIndex, buf, BUFSIZE, curPos, pos);
		line = nextLine(in, gzIndex, buf, BUFSIZE, curPos, pos);
		line = nextLine(in, gzIndex, buf, BUFSIZE, curPos, pos);

		oldPos = pos;
		line = nextLine(in, gzIndex, buf, BUFSIZE, curPos, pos);
	}

	free(buf);

	return barcodesIndex;
}



BarcodesIndex indexBarcodesFromFastqGz(string fastqFile, unsigned nbThreads) {
	// Open the necessary number of streams
	vector<FILE*> ins(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		ins[i] = fopen(fastqFile.c_str(), "rb");
		if (ins[i] == NULL) {
			throw ios_base::failure("open: Unable to open gziped fastq file " + fastqFile + ". Please provide an existing and valid file.");
		}
	}

	struct access* gzIndex = (struct access*) malloc(sizeof(struct access));
	gzIndex = deserializeGzIndex(gzIndex, fastqFile + "i");

	// Split the processing of the fastq file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<BarcodesIndex>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(indexBarcodesFromFastqGzOffsets, ins[i], gzIndex, i * gzIndex->size / nbThreads, min((int) ((i + 1) * gzIndex->size / nbThreads), gzIndex->size - 1), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	BarcodesIndex barcodesIndex;
	BarcodesIndex curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			for (auto v : p.second) {
				barcodesIndex[p.first].push_back(v);
			}
		}
	}

	// Close streams
	for (unsigned i = 0; i < nbThreads; i++) {
		fclose(ins[i]);
	}

	// Free index
	freeGzIndex(gzIndex);

	return barcodesIndex;
}

void saveBarcodesIndex(BarcodesIndex& barcodesIndex, string file) {
	ofstream out;
	out.open(file, ios::out | ios::binary);
	if (!out.is_open()) {
		throw ios_base::failure("open: Unable to open file " + file + ".");
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
		throw ios_base::failure("open: Unable to open barcodes index file " + file + ". Please provide an existing and valid file.");
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