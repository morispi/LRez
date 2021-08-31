#ifndef __LREZ_INDEX_MANAGEMENT_FASTQ__
#define __LREZ_INDEX_MANAGEMENT_FASTQ__

#include <vector>
#include "robin_hood.h"
#include "utils.h"

using namespace std;

#define BUFSIZE 134217728

typedef robin_hood::unordered_map<barcode, vector<int64_t>> BarcodesIndex;


/**
	Index the offsets of the barcodes contained in a given fastq file.

	@param fastqFile fastq file to build the index from
	@param nbThreads number of threads to use, set to 1 by default
	@throws ios_base::failure thrown if fastqFile could not be open
	@return the barcode index associating barcodes to the set of offsets they occur at
*/
BarcodesIndex indexBarcodesFromFastq(string fastqFile, unsigned nbThreads = 1);


/**
	Index the offsets of the barcodes contained in a given gzipped fastq file.

	@param fastqFile gzipped fastq file to build the index from
	@param nbThreads number of threads to use, set to 1 by default
	@throws ios_base::failure thrown if fastqFile could not be open
	@return the barcode index associating barcodes to the set of offsets they occur at
*/
BarcodesIndex indexBarcodesFromFastqGz(string fastqFile, unsigned nbThreads = 1);

/**
	Save the content of the barcode index in a given file.

	@param index barcode offsets index to save
	@throws ios_base::failure thrown if file could not be open
	@param file file where to store the index
*/
void saveBarcodesIndex(BarcodesIndex& index, string file);

/**
	Load the barcode index stored in a given file.

	@param file file where the barcode index is stored
	@throws ios_base::failure thrown if file could not be open
	@return the barcode index filled with the data contained in the file
*/
BarcodesIndex loadBarcodesIndex(string file);

#endif