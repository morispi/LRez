#include "alignmentsRetrieval.h"
#include "reverseComplement.h"
#include "../CTPL/ctpl_stl.h"

vector<BamAlignment> retrieveAlignmentsWithBarcodeBits_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, barcode b) {
	vector<BamAlignment> res;
	BamAlignment al;

	for (int64_t r : BarcodesOffsetsIndex[b]) {
		if (!reader.FSeek(r)) {
			throw runtime_error("Fseek: Error while attempting to jump to offset " + to_string(r) + ".");
		}
		
		reader.GetNextAlignment(al);
		res.push_back(al);
	}

	return res;
}

vector<BamAlignment> retrieveAlignmentsWithBarcode_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string bc) {
	return retrieveAlignmentsWithBarcodeBits_BamReader(reader, BarcodesOffsetsIndex, stringToBarcode(bc));
}

vector<BamAlignment> retrieveAlignmentsWithBarcodeBits(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, barcode bc) {
	BamReader reader;
	if (!reader.Open(bamFile)) {
		throw ios_base::failure("Open: Unable to open BAM file " + bamFile + ". Please make sure the file exists.");
	}
	if (!reader.LocateIndex()) {
		throw ios_base::failure("LocateIndex: Unable to find a BAM index for file " + bamFile + ". Please build the BAM index or provide a BAM file for which the BAM index is built.");
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
		throw ios_base::failure("open: Unable to open barcodes list file " + barcodesList + ". Please provide an existing and valid file.");
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

/**
	Querie the index to retrieve alignments that have their barcodes in the specified list.
*/
vector<BamAlignment> retrieveAlignmentsWithBarcodesList(int id, BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, vector<string> queries) {
	vector<BamAlignment> res;
	vector<BamAlignment> tmpRes;
	for (string q : queries) {
		tmpRes = retrieveAlignmentsWithBarcode_BamReader(reader, BarcodesOffsetsIndex, q);
		res.insert(res.end(), tmpRes.begin(), tmpRes.end());
	}

	return res;
}

vector<BamAlignment> retrieveAlignmentsWithBarcodes(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string barcodesList, unsigned nbThreads) {
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
	vector<std::future<vector<BamAlignment>>> results(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		results[i] = myPool.push(retrieveAlignmentsWithBarcodesList, ref(readers[i]), ref(BarcodesOffsetsIndex), ref(queries[i]));
	}

	// Retrieve threads subresults and build global result
	vector<BamAlignment> res;
	vector<BamAlignment> curRes;
	for (unsigned i = 0; i < nbThreads; i++) {
		curRes = results[i].get();
		res.insert(res.end(), curRes.begin(), curRes.end());
	}

	// Close BamReaders
	for (unsigned i = 0; i < nbThreads; i++) {
		readers[i].Close();
	}

	return res;

}