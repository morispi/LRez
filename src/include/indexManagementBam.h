#ifndef __LREZ_INDEX_MANAGEMENT_BAM__
#define __LREZ_INDEX_MANAGEMENT_BAM__

#include <vector>
#include "robin_hood.h"
#include "api/BamAux.h"
#include "utils.h"

using namespace std;
using namespace BamTools;

typedef robin_hood::unordered_map<barcode, vector<int64_t>> BarcodesOffsetsIndex;

typedef robin_hood::unordered_map<barcode, vector<pair<int32_t, int32_t>>> BarcodesPositionsIndex;


/**
	Index the (chromosome, begPosition) pairs of the barcodes contained in a given BAM file.

	@param reader bamFile BAM file to build the index from
	@param primary whether to only index primary alignments or not
	@param quality minimum quality to index an alignment
	@return the barcode positions index associating barcodes to the set of (chromosome, begPosition) pairs they occur at
*/
BarcodesPositionsIndex indexBarcodesPositionsFromBam(string bamFile, bool primary, unsigned quality);

/**
	Save the content of the barcode positions index in a given file.

	@param BarcodesPositionsIndex barcode positions index to save
	@param file file where to store the index
*/
void saveBarcodesPositionsIndex(BarcodesPositionsIndex& index, string file);

/**
	Load the barcode positions index stored in a given file.

	@param file file where the barcode positions index is stored
	@return the barcode positions index filled with the data contained in the file
*/
BarcodesPositionsIndex loadBarcodesPositionsIndex(string file);


/**
	Index the offsets of the barcodes contained in a given BAM file.

	@param reader bamFile BAM file to build the index from
	@param primary whether to only index primary alignments or not
	@param quality minimum quality to index an alignment
	@return the barcode index associating barcodes to the set of offsets they occur at
*/
BarcodesOffsetsIndex indexBarcodesOffsetsFromBam(string bamFile, bool primary, unsigned quality);

/**
	Save the content of the barcode index in a given file.

	@param BarcodesOffsetsIndex barcode offsets index to save
	@param file file where to store the index
*/
void saveBarcodesOffsetsIndex(BarcodesOffsetsIndex& index, string file);

/**
	Load the barcode index stored in a given file.

	@param file file where the barcode offsets index is stored
	@return the barcode offsets index filled with the data contained in the file
*/
BarcodesOffsetsIndex loadBarcodesOffsetsIndex(string file);

#endif