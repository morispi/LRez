#include "alignmentsRetrieval.h"
#include "reverseComplement.h"
#include <chrono>


using namespace std::chrono;

BamAlignment retrieveAlignmentWithBarcode(BamReader& reader, int64_t position, barcode b) {
	BamAlignment al;
	string tag;

	reader.GetNextAlignment(al);
	return al;
}

vector<BamAlignment> retrieveAlignmentsWithBarcodeBits_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, barcode b) {
	vector<BamAlignment> res;
	BamAlignment al;

	for (int64_t r : BarcodesOffsetsIndex[b]) {
		if (!reader.FSeek(r)) {
			fprintf(stderr, "Error while attempting to jump to offset %ld", r);
			exit(EXIT_FAILURE);
		}
		
		res.push_back(retrieveAlignmentWithBarcode(reader, r, b));
	}

	return res;
}

vector<BamAlignment> retrieveAlignmentsWithBarcode_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string bc) {
	return retrieveAlignmentsWithBarcodeBits_BamReader(reader, BarcodesOffsetsIndex, stringToBarcode(bc));
}

vector<BamAlignment> retrieveAlignmentsWithBarcodeBits(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, barcode bc) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	
	vector<BamAlignment> res = retrieveAlignmentsWithBarcodeBits_BamReader(reader, BarcodesOffsetsIndex, bc);
	reader.Close();

	return res;
}

vector<BamAlignment> retrieveAlignmentsWithBarcode(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string bc) {
	return retrieveAlignmentsWithBarcodeBits(bamFile, BarcodesOffsetsIndex, stringToBarcode(bc));
}

vector<BamAlignment> retrieveAlignmentsWithBarcodes_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string barcodesList) {
	vector<BamAlignment> res, tmp;
	string line;	

	ifstream bc;
	bc.open(barcodesList);
	if (!bc.is_open()) {
		fprintf(stderr, "Unable to open barcodes list file %s. Please provide an existing and valid file.\n", barcodesList.c_str());
		exit(EXIT_FAILURE);
	}

	while (getline(bc, line)) {
		tmp = retrieveAlignmentsWithBarcode_BamReader(reader, BarcodesOffsetsIndex, line);
		for (BamAlignment al : tmp) {
			res.push_back(al);
		}
	}

	bc.close();

	return res;
}

vector<BamAlignment> retrieveAlignmentsWithBarcodes(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string barcodesList) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	
	vector<BamAlignment> res = retrieveAlignmentsWithBarcodes_BamReader(reader, BarcodesOffsetsIndex, barcodesList);
	reader.Close();

	return res;
}