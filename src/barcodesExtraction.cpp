#include "barcodesExtraction.h"
#include "reverseComplement.h"
#include "../CTPL/ctpl_stl.h"

robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		throw runtime_error("stringToBamRegion: Unable to parse region " + region + ". Please ensure the region of interest is in format chromosome:startPosition-endPosition.");
	}
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + region + ".");
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
		throw runtime_error("stringToBamRegion: Unable to parse region " + region + ". Please ensure the region of interest is in format chromosome:startPosition-endPosition.");
	}
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + region + ".");
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
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
	}

	robin_hood::unordered_set<barcode> res = extractBarcodesBitsFromRegion_BamReader(reader, region);
	reader.Close();

	return res;

}

robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion(string bamFile, string region) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
	}

	robin_hood::unordered_set<string> res = extractBarcodesSeqsFromRegion_BamReader(reader, region);
	reader.Close();

	return res;
}

/**
	Extract the barcodes (binary representation) from the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromBamRegions(int id, BamReader& reader, string begRegion, string endRegion, bool tilEnd = 0) {
	robin_hood::unordered_set<barcode> barcodes;
	BamAlignment al;
	barcode b;
	string bc, barcode;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + begRegion + " ; " + endRegion + ".");
	}

	// Process and extract barcodes from alignments located between the two regions of interest
	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodes.insert(b);
			}
		}
		pos = reader.FTell();
	}

	// Treat unmapped reads if we're at the end of the file
	if (tilEnd) {
		reader.Rewind();
		if (!reader.FSeek(pos)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(pos) + ".");
		}

		while (reader.GetNextAlignment(al)) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodes.insert(b);
			}
			pos = reader.FTell();
		}
	}

	return barcodes;
}

robin_hood::unordered_set<barcode> extractBarcodesBitsFromBAM(string bamFile, unsigned nbThreads) {
	// Open the necessary number of BamReaders
	vector<BamReader> readers(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		if (!readers[i].Open(bamFile)) {
			throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
		}
		if (!readers[i].LocateIndex()) {
			throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
		}
		readers[i].Rewind();
	}

	// Split the reference genome into regions of size 10000
	vector<string> regions = extractRegionsList(readers[0], 10000);
	readers[0].Rewind();

	// Split the processing of the BAM file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<robin_hood::unordered_set<barcode>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(extractBarcodesBitsFromBamRegions, ref(readers[i]), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	robin_hood::unordered_set<barcode> barcodes;
	robin_hood::unordered_set<barcode> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			barcodes.insert(p);
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return barcodes;
}

/**
	Extract the barcodes (string representation) from the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromBamRegions(int id, BamReader& reader, string begRegion, string endRegion, bool tilEnd = 0) {
	robin_hood::unordered_set<string> barcodes;
	BamAlignment al;
	string bc, barcode;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + begRegion + " ; " + endRegion + ".");
	}

	// Process and extract barcodes from alignments located between the two regions of interest
	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.insert(barcode);
			}
		}
		pos = reader.FTell();
	}

	// Treat unmapped reads if we're at the end of the file
	if (tilEnd) {
		reader.Rewind();
		if (!reader.FSeek(pos)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(pos) + ".");
		}
		while (reader.GetNextAlignment(al)) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.insert(barcode);
			}
			pos = reader.FTell();
		}
	}

	return barcodes;
}

robin_hood::unordered_set<string> extractBarcodesSeqsFromBAM(string bamFile, unsigned nbThreads) {
	// Open the necessary number of BamReaders
	vector<BamReader> readers(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		if (!readers[i].Open(bamFile)) {
			throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
		}
		if (!readers[i].LocateIndex()) {
			throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
		}
		readers[i].Rewind();
	}

	// Split the reference genome into regions of size 10000
	vector<string> regions = extractRegionsList(readers[0], 10000);
	readers[0].Rewind();

	// Split the processing of the BAM file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<robin_hood::unordered_set<string>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(extractBarcodesSeqsFromBamRegions, ref(readers[i]), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	robin_hood::unordered_set<string> barcodes;
	robin_hood::unordered_set<string> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			barcodes.insert(p);
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return barcodes;
}

vector<barcode> extractBarcodesBitsFromRegionWithDuplicates_BamReader(BamReader& reader, string region) {
	BamRegion r = stringToBamRegion(reader, region);
	if (r.isNull()) {
		throw runtime_error("stringToBamRegion: Unable to parse region " + region + ". Please ensure the region of interest is in format chromosome:startPosition-endPosition.");
	}
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + region + ".");
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
		throw runtime_error("stringToBamRegion: Unable to parse region " + region + ". Please ensure the region of interest is in format chromosome:startPosition-endPosition.");
	}
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + region + ".");
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
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
	}

	vector<barcode> res = extractBarcodesBitsFromRegionWithDuplicates_BamReader(reader, region);
	reader.Close();

	return res;

}

vector<string> extractBarcodesSeqsFromRegionWithDuplicates(string bamFile, string region) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
	}

	vector<string> res = extractBarcodesSeqsFromRegionWithDuplicates_BamReader(reader, region);
	reader.Close();

	return res;
}

/**
	Extract the barcodes (binary representation) from the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
vector<barcode> extractBarcodesBitsFromBamRegionsWithDuplicates(int id, BamReader& reader, string begRegion, string endRegion, bool tilEnd = 0) {
	vector<barcode> barcodes;
	BamAlignment al;
	barcode b;
	string bc, barcode;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + begRegion + " ; " + endRegion + ".");
	}

	// Process and extract barcodes from alignments located between the two regions of interest
	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodes.push_back(b);
			}
		}
		pos = reader.FTell();
	}

	// Treat unmapped reads if we're at the end of the file
	if (tilEnd) {
		reader.Rewind();
		if (!reader.FSeek(pos)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(pos) + ".");
		}

		while (reader.GetNextAlignment(al)) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				b = stringToBarcode(barcode);
				barcodes.push_back(b);
			}
			pos = reader.FTell();
		}
	}

	return barcodes;
}

vector<barcode> extractBarcodesBitsFromBAMWithDuplicates(string bamFile, unsigned nbThreads) {
	// Open the necessary number of BamReaders
	vector<BamReader> readers(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		if (!readers[i].Open(bamFile)) {
			throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
		}
		if (!readers[i].LocateIndex()) {
			throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
		}
		readers[i].Rewind();
	}

	// Split the reference genome into regions of size 10000
	vector<string> regions = extractRegionsList(readers[0], 10000);
	readers[0].Rewind();

	// Split the processing of the BAM file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<vector<barcode>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(extractBarcodesBitsFromBamRegionsWithDuplicates, ref(readers[i]), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	vector<barcode> barcodes;
	vector<barcode> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			barcodes.push_back(p);
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return barcodes;
}

/**
	Extract the barcodes (string representation) from the portion of the BAM file located between begRegion and endRegion.
	tilEnd has to be set to one for the thread processing the end of the BAM file, so it can process unmapped reads.
*/
vector<string> extractBarcodesSeqsFromBamRegionsWithDuplicates(int id, BamReader& reader, string begRegion, string endRegion, bool tilEnd = 0) {
	vector<string> barcodes;
	BamAlignment al;
	string bc, barcode;

	// Transform the two regions into a single bam region, and jump to this region of the BAM file
	vector<string> t1 = splitString(begRegion, ":");
	vector<string> tt1 = splitString(t1[1], "-");
	vector<string> t2 = splitString(endRegion, ":");
	vector<string> tt2 = splitString(t2[1], "-");
	BamRegion r(reader.GetReferenceID(t1[0]), stoi(tt1[0]), reader.GetReferenceID(t2[0]), stoi(tt2[1]));
	if (!reader.SetRegion(r)) {
		throw runtime_error("Error while attempting to jump to region " + begRegion + " ; " + endRegion + ".");
	}

	// Process and extract barcodes from alignments located between the two regions of interest
	int64_t pos = reader.FTell();
	while (reader.GetNextAlignment(al)) {
		bool processAlignment = (al.RefID > reader.GetReferenceID(t1[0]) or (al.RefID == reader.GetReferenceID(t1[0]) and al.Position >= stoi(tt1[0])))
		and ( (al.RefID < reader.GetReferenceID(t2[0]) or (al.RefID == reader.GetReferenceID(t2[0]) and (al.Position < stoi(tt2[0]) or tilEnd))));
		if (processAlignment) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.push_back(barcode);
			}
		}
		pos = reader.FTell();
	}

	// Treat unmapped reads if we're at the end of the file
	if (tilEnd) {
		reader.Rewind();
		if (!reader.FSeek(pos)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(pos) + ".");
		}
		while (reader.GetNextAlignment(al)) {
			barcode.clear();
			al.GetTag(BXTAG, barcode);
			if (isValidBarcode(barcode)) {
				barcodes.push_back(barcode);
			}
			pos = reader.FTell();
		}
	}

	return barcodes;
}

vector<string> extractBarcodesSeqsFromBAMWithDuplicates(string bamFile, unsigned nbThreads) {
	// Open the necessary number of BamReaders
	vector<BamReader> readers(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		if (!readers[i].Open(bamFile)) {
			throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
		}
		if (!readers[i].LocateIndex()) {
			throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
		}
		readers[i].Rewind();
	}

	// Split the reference genome into regions of size 10000
	vector<string> regions = extractRegionsList(readers[0], 10000);
	readers[0].Rewind();

	// Split the processing of the BAM file into separate threads
	ctpl::thread_pool myPool(nbThreads);
	vector<std::future<vector<string>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(extractBarcodesSeqsFromBamRegionsWithDuplicates, ref(readers[i]), ref(regions[i * regions.size() / nbThreads]), ref(regions[min((i + 1) * regions.size() / nbThreads, regions.size() - 1)]), i == nbThreads - 1);
	}

	// Retrieve threads subresults and build global index
	vector<string> barcodes;
	vector<string> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		for (auto p : curRes) {
			barcodes.push_back(p);
		}
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return barcodes;
}