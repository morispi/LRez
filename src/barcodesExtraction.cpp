#include "barcodesExtraction.h"
#include "reverseComplement.h"

robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		fprintf(stderr, "Unable to parse region %s. Please ensure the region of interest is in format chromosome:startPosition-endPosition.\n", region.c_str());
		exit(EXIT_FAILURE);
	}
	if (!reader.SetRegion(r)) {
		fprintf(stderr, "Error while attempting to jump to region %d:%d-%d:%d", r.LeftRefID, r.LeftPosition, r.RightRefID, r.RightPosition);
		exit(EXIT_FAILURE);
	}
	

	robin_hood::unordered_set<barcode> barcodes;
	BamAlignment al;
	string bc, barcode;
	while (reader.GetNextAlignment(al) and al.RefID <= r.RightRefID and al.Position <= r.RightPosition) {
		if (al.IsMapped()) {
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
				barcodes.insert(stringToBarcode(barcode));
			}
		}
	}

	reader.Rewind();
	return barcodes;
}

robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion_BamReader(BamReader& reader, string region) {
	robin_hood::unordered_set<barcode> bc = extractBarcodesBitsFromRegion_BamReader(reader, region);
	
	robin_hood::unordered_set<string> barcodes;
	for (auto b : bc) {
		barcodes.insert(barcodeToString(b));
	}

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
	robin_hood::unordered_set<barcode> bc = extractBarcodesBitsFromRegion(bamFile, region);
	
	robin_hood::unordered_set<string> barcodes;
	for (auto b : bc) {
		barcodes.insert(barcodeToString(b));
	}

	return barcodes;
}

robin_hood::unordered_set<barcode> extractBarcodesBitsFromBAM(string bamFile) {
	// if (CONSIDER_RX) {
	// 	cerr << "considering RX" << endl;
	// } else {
	// 	cerr << "nope" << endl;
	// }
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
			barcodes.insert(stringToBarcode(barcode));
		}
	}

	reader.Close();
	return barcodes;
}


robin_hood::unordered_set<string> extractBarcodesSeqsFromBAM(string bamFile) {
	robin_hood::unordered_set<barcode> bc = extractBarcodesBitsFromBAM(bamFile);

	robin_hood::unordered_set<string> barcodes;
	for (auto b : bc) {
		barcodes.insert(barcodeToString(b));
	}

	return barcodes;
}