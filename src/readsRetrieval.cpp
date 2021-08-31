#include "readsRetrieval.h"
#include "reverseComplement.h"
#include "../CTPL/ctpl_stl.h"

vector<string> retrieveReadsWithBarcodeBits_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, barcode bc) {
	vector<string> res;
	string line, read;

	for (int64_t r : BarcodesIndex[bc]) {
		if (!in.seekg(r, in.beg)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(r) + ".");
		}
		
		getline(in, line);
		read = line;
		getline(in, line);
		read = read + "\n" + line;
		getline(in, line);
		read = read + "\n" + line;
		getline(in, line);
		read = read + "\n" + line;
		
		res.push_back(read);
	}

	return res;
}

vector<string> retrieveReadsWithBarcodeBits(string fastqFile, BarcodesIndex& BarcodesIndex, barcode bc) {
	ifstream in;
	in.open(fastqFile);
	if (!in.is_open()) {
		throw ios_base::failure("open: Unable to open fastq file " + fastqFile + ". Please provide an existing and valid file.");
	}
	
	vector<string> res = retrieveReadsWithBarcodeBits_Stream(in, BarcodesIndex, bc);
	in.close();

	return res;
}


vector<string> retrieveReadsWithBarcode_Stream(ifstream& reader, BarcodesIndex& BarcodesIndex, string bc) {
	return retrieveReadsWithBarcodeBits_Stream(reader, BarcodesIndex, stringToBarcode(bc));
}

vector<string> retrieveReadsWithBarcode(string fastqFile, BarcodesIndex& BarcodesIndex, string bc) {
	return retrieveReadsWithBarcodeBits(fastqFile, BarcodesIndex, stringToBarcode(bc));
}

vector<string> retrieveReadsWithBarcodes_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, string barcodesList) {
	vector <string> res, tmp;
	string line;

	ifstream bc;
	bc.open(barcodesList);
	if (!bc.is_open()) {
		throw ios_base::failure("open: Unable to open barcodes list file " + barcodesList + ". Please provide an existing and valid file.");
	}

	while (getline(bc, line)) {
		tmp = retrieveReadsWithBarcode_Stream(in, BarcodesIndex, line);
		for (string s : tmp) {
			res.push_back(s);
		}
	}

	bc.close();

	return res;
}

/**
	Query the index to retrieve alignments that have their barcodes in the specified list.
*/
vector<string> retrieveReadsWithBarcodesList(int id, ifstream& in, BarcodesIndex& BarcodesIndex, vector<string> queries) {
	vector<string> res;
	vector<string> tmpRes;
	for (string q : queries) {
		tmpRes = retrieveReadsWithBarcode_Stream(in, BarcodesIndex, q);
		res.insert(res.end(), tmpRes.begin(), tmpRes.end());
	}

	return res;
}

vector<string> retrieveReadsWithBarcodes(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList, unsigned nbThreads) {
	// Open the necessary number of streams
	vector<ifstream> ins(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		ins[i].open(fastqFile);
		if (!ins[i].is_open()) {
			throw ios_base::failure("open: Unable to open fastq file " + fastqFile + ". Please provide an existing and valid file.");
		}
	}

	// Fill vectors with the query barcodes
	ifstream bc;
	bc.open(barcodesList);
	if (!bc.is_open()) {
		throw ios_base::failure("open: Unable to open barcodes list file " + barcodesList + ". Please provide an existing and valid file.");
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
	vector<std::future<vector<string>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(retrieveReadsWithBarcodesList, ref(ins[i]), ref(BarcodesIndex), ref(queries[i]));
	}

	// Retrieve threads subresults and build global result
	vector<string> res;
	vector<string> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		res.insert(res.end(), curRes.begin(), curRes.end());
	}

	// Close streams
	for (unsigned i = 0; i < nbThreads; i++) {
		ins[i].close();
	}

	return res;
}

vector<string> retrieveReadsWithBarcodeBits_Gzip_Stream_Index(FILE* in, struct access* gzIndex, BarcodesIndex& barcodesIndex, barcode bc) {
	vector<string> res;
	string read;

	for (int64_t r : barcodesIndex[bc]) {
		read = extractFastqReadFromOffset(in, gzIndex, r);
		res.push_back(read);
	}

	return res;
}

vector<string> retrieveReadsWithBarcodeBits_Gzip(string fastqFile, BarcodesIndex& barcodesIndex, barcode bc) {
	FILE *in;
	in = fopen(fastqFile.c_str(), "rb");
	if (in == NULL) {
		throw ios_base::failure("open: Unable to open gziped fastq file " + fastqFile + ". Please provide an existing and valid file.");
	}

	struct access* gzIndex = NULL;
	gzIndex = deserializeGzIndex(gzIndex, fastqFile + "i");
	
	vector<string> res = retrieveReadsWithBarcodeBits_Gzip_Stream_Index(in, gzIndex, barcodesIndex, bc);
	fclose(in);
	freeGzIndex(gzIndex);

	return res;
}

vector<string> retrieveReadsWithBarcode_Gzip_Stream_Index(FILE* in, struct access* gzIndex, BarcodesIndex& BarcodesIndex, string barcode) {
	return retrieveReadsWithBarcodeBits_Gzip_Stream_Index(in, gzIndex, BarcodesIndex, stringToBarcode(barcode));
}

vector<string> retrieveReadsWithBarcode_Gzip(string fastqFile, BarcodesIndex& BarcodesIndex, string barcode) {
	return retrieveReadsWithBarcodeBits_Gzip(fastqFile, BarcodesIndex, stringToBarcode(barcode));
}

vector<string> retrieveReadsWithBarcodes_Gzip_Stream_Index(FILE* in, struct access* gzIndex, BarcodesIndex& BarcodesIndex, string barcodesList) {
	vector <string> res, tmp;
	string line;

	ifstream bc;
	bc.open(barcodesList);
	if (!bc.is_open()) {
		throw ios_base::failure("open: Unable to open barcodes list file " + barcodesList + ". Please provide an existing and valid file.");
	}

	while (getline(bc, line)) {
		tmp = retrieveReadsWithBarcode_Gzip_Stream_Index(in, gzIndex, BarcodesIndex, line);
		for (string s : tmp) {
			res.push_back(s);
		}
	}

	bc.close();

	return res;
}

/**
	Querie the index to retrieve alignments that have their barcodes in the specified list.
*/
vector<string> retrieveReadsWithBarcodesList_Gzip(int id, FILE* in, struct access* gzIndex, BarcodesIndex& BarcodesIndex, vector<string> queries) {
	vector<string> res;
	vector<string> tmpRes;
	for (string q : queries) {
		tmpRes = retrieveReadsWithBarcode_Gzip_Stream_Index(in, gzIndex, BarcodesIndex, q);
		res.insert(res.end(), tmpRes.begin(), tmpRes.end());
	}

	return res;
}

vector<string> retrieveReadsWithBarcodes_Gzip(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList, unsigned nbThreads) {
	// Open the necessary number of streams
	vector<FILE*> ins(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		ins[i] = fopen(fastqFile.c_str(), "rb");
		if (ins[i] == NULL) {
			throw ios_base::failure("open: Unable to open gziped fastq file " + fastqFile + ". Please provide an existing and valid file.");
		}
	}

	vector<struct access*> gzIndexes(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		gzIndexes[i] = deserializeGzIndex(gzIndexes[i], fastqFile + "i");
	}

	// Fill vectors with the query barcodes
	ifstream bc;
	bc.open(barcodesList);
	if (!bc.is_open()) {
		throw ios_base::failure("open: Unable to open barcodes list file " + barcodesList + ". Please provide an existing and valid file.");
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
	vector<std::future<vector<string>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(retrieveReadsWithBarcodesList_Gzip, ins[i], gzIndexes[i], ref(BarcodesIndex), ref(queries[i]));
	}

	// Retrieve threads subresults and build global result
	vector<string> res;
	vector<string> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		res.insert(res.end(), curRes.begin(), curRes.end());
	}

	// Close streams
	for (unsigned i = 0; i < nbThreads; i++) {
		fclose(ins[i]);
	}

	// Free indexes
	for (unsigned i = 0; i < nbThreads; i++) {
		freeGzIndex(gzIndexes[i]);
	}

	return res;
}