#include "barcodesExtraction.h"
#include "reverseComplement.h"

robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		fprintf(stderr, "Unable to parse region %s. Please ensure the region of interest is in format chromosome:startPosition-endPosition.\n", region.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.SetRegion(r)) {
		fprintf(stderr, "Error while attempting to jump to region %s", region.c_str());
		exit(EXIT_FAILURE);
	}
	

	robin_hood::unordered_set<barcode> barcodes;
	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al) and al.RefID <= r.RightRefID and al.Position <= r.RightPosition) {
		if (al.IsMapped()) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.insert(stringToBarcode(barcode));
			}
		}
	}

	reader.Rewind();
	return barcodes;
}

robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		fprintf(stderr, "Unable to parse region %s. Please ensure the region of interest is in format chromosome:startPosition-endPosition.\n", region.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.SetRegion(r)) {
		fprintf(stderr, "Error while attempting to jump to region %s", region.c_str());
		exit(EXIT_FAILURE);
	}
	

	robin_hood::unordered_set<string> barcodes;
	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al) and al.RefID <= r.RightRefID and al.Position <= r.RightPosition) {
		if (al.IsMapped()) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.insert(barcode);
			}
		}
	}

	reader.Rewind();
	return barcodes;
}

robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion(string bamFile, string region) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	robin_hood::unordered_set<barcode> res = extractBarcodesBitsFromRegion_BamReader(reader, region);
	reader.Close();

	return res;

}

robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion(string bamFile, string region) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	robin_hood::unordered_set<string> res = extractBarcodesSeqsFromRegion_BamReader(reader, region);
	reader.Close();

	return res;
}

robin_hood::unordered_set<barcode> extractBarcodesBitsFromBAM(string bamFile) {
	robin_hood::unordered_set<barcode> barcodes;
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al)) {
		barcode.clear();
		al.GetTag(BXTAG, barcode);
		if (isValidBarcode(barcode)) {
			barcodes.insert(stringToBarcode(barcode));
		}
	}

	reader.Close();
	return barcodes;
}


robin_hood::unordered_set<string> extractBarcodesSeqsFromBAM(string bamFile) {
	robin_hood::unordered_set<string> barcodes;
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al)) {
		barcode.clear();
		al.GetTag(BXTAG, barcode);
		if (isValidBarcode(barcode)) {
			barcodes.insert(barcode);
		}
	}

	reader.Close();
	return barcodes;
}

vector<barcode> extractBarcodesBitsFromRegionWithDuplicates_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		fprintf(stderr, "Unable to parse region %s. Please ensure the region of interest is in format chromosome:startPosition-endPosition.\n", region.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.SetRegion(r)) {
		fprintf(stderr, "Error while attempting to jump to region %s", region.c_str());
		exit(EXIT_FAILURE);
	}
	

	vector<barcode> barcodes;
	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al) and al.RefID <= r.RightRefID and al.Position <= r.RightPosition) {
		if (al.IsMapped()) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.push_back(stringToBarcode(barcode));
			}
		}
	}

	reader.Rewind();
	return barcodes;
}

vector<string> extractBarcodesSeqsFromRegionWithDuplicates_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		fprintf(stderr, "Unable to parse region %s. Please ensure the region of interest is in format chromosome:startPosition-endPosition.\n", region.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.SetRegion(r)) {
		fprintf(stderr, "Error while attempting to jump to region %s", region.c_str());
		exit(EXIT_FAILURE);
	}
	

	vector<string> barcodes;
	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al) and al.RefID <= r.RightRefID and al.Position <= r.RightPosition) {
		if (al.IsMapped()) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.push_back(barcode);
			}
		}
	}

	reader.Rewind();
	return barcodes;
}

vector<barcode> extractBarcodesBitsFromRegionWithDuplicates(string bamFile, string region) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	vector<barcode> res = extractBarcodesBitsFromRegionWithDuplicates_BamReader(reader, region);
	reader.Close();

	return res;

}

vector<string> extractBarcodesSeqsFromRegionWithDuplicates(string bamFile, string region) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.LocateIndex()) {
		fprintf(stderr, "Unable to find a BAM index for file %s. Please build the BAM index or provide a BAM file for which the BAM index is built\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	vector<string> res = extractBarcodesSeqsFromRegionWithDuplicates_BamReader(reader, region);
	reader.Close();

	return res;
}

vector<barcode> extractBarcodesBitsFromBAMWithDuplicates(string bamFile) {
	vector<barcode> barcodes;
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al)) {
		barcode.clear();
		al.GetTag(BXTAG, barcode);
		if (isValidBarcode(barcode)) {
			barcodes.push_back(stringToBarcode(barcode));
		}
	}

	reader.Close();
	return barcodes;
}


vector<string> extractBarcodesSeqsFromBAMWithDuplicates(string bamFile) {
	vector<string> barcodes;
	BamReader reader;
	if (!reader.Open(bamFile)) {
		fprintf(stderr, "Unable open BAM file %s. Please make sure the file exists.\n", bamFile.c_str());
		exit(EXIT_FAILURE);
	}

	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al)) {
		barcode.clear();
		al.GetTag(BXTAG, barcode);
		if (isValidBarcode(barcode)) {
			barcodes.push_back(barcode);
		}
	}

	reader.Close();
	return barcodes;
}