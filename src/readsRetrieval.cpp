#include "readsRetrieval.h"
#include "reverseComplement.h"

vector<string> retrieveReadsWithBarcodeBits_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, barcode bc) {
	vector<string> res;
	string line, read;

	for (int64_t r : BarcodesIndex[bc]) {
		if (!in.seekg(r, in.beg)) {
			fprintf(stderr, "Error while attempting to jump to offset %ld\n", r);
			exit(EXIT_FAILURE);
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
		fprintf(stderr, "Unable to open barcodes index file %s. Please provide an existing and valid file.\n", fastqFile.c_str());
		exit(EXIT_FAILURE);
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
		fprintf(stderr, "Unable to open barcodes list file %s. Please provide an existing and valid file.\n", barcodesList.c_str());
		exit(EXIT_FAILURE);
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

vector<string> retrieveReadsWithBarcodes(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList) {
	ifstream in;
	in.open(fastqFile);
	if (!in.is_open()) {
		fprintf(stderr, "Unable to open fastq file %s. Please provide an existing and valid file.\n", fastqFile.c_str());
		exit(EXIT_FAILURE);
	}

	vector<string> res = retrieveReadsWithBarcodes_Stream(in, BarcodesIndex, barcodesList);
	in.close();

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
		fprintf(stderr, "Unable to open gzip fastq file %s. Please provide an existing and valid file.\n", fastqFile.c_str());
	  	exit(EXIT_FAILURE);
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
		fprintf(stderr, "Unable to open barcodes list file %s. Please provide an existing and valid file.\n", barcodesList.c_str());
		exit(EXIT_FAILURE);
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

vector<string> retrieveReadsWithBarcodes_Gzip(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList) {
	FILE *in;
	in = fopen(fastqFile.c_str(), "rb");
	if (in == NULL) {
		fprintf(stderr, "Unable to open gzip fastq file %s. Please provide an existing and valid file.\n", fastqFile.c_str());
	  	exit(EXIT_FAILURE);
	}

	struct access* gzIndex = NULL;
	gzIndex = deserializeGzIndex(gzIndex, fastqFile + "i");

	vector<string> res = retrieveReadsWithBarcodes_Gzip_Stream_Index(in, gzIndex, BarcodesIndex, barcodesList);
	fclose(in);
	freeGzIndex(gzIndex);

	return res;
}