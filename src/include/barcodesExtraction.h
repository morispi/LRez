#ifndef __LREZ_BARCODES_EXTRACTION__
#define __LREZ_BARCODES_EXTRACTION__

#include "utils.h"
#include "robin_hood.h"

/**
	Extract the barcodes of a given region, in binary representation.

	@param reader BamReader open on the desired BAM file
	@param region region of interest
	@return a vector containing the barcodes present in the region
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion_BamReader(BamReader& reader, string region);

/**
	Extract the barcodes of a given region, in string representation.

	@param reader BamReader open on the desired BAM file
	@param region region of interest
	@return a vector containing the barcodes present in the region
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion_BamReader(BamReader& reader, string region);

/**
	Extract the barcodes of a given region, in binary representation.

	@param reader bamFile BAM file to extract barcodes from
	@param region region of interest
	@return a vector containing the barcodes present in the region, in binary representation
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion(string bamFile, string region);


/**
	Extract the barcodes of a given region, in string representation.

	@param reader bamFile BAM file to extract barcodes from
	@param region region of interest
	@return a vector containing the barcodes present in the region, in binary string
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion(string bamFile, string region);

/**
	Extract all the barcodes of a given BAM file, in binary representation.

	@param bamFile bamFile BAM file to extract barcodes from
	@return a vector containing the barcodes present in the BAM file, in binary representation
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromBAM(string bamFile);

/**
	Extract all the barcodes of a given BAM file, in string representation.

	@param bamFile bamFile BAM file to extract barcodes from
	@return a vector containing the barcodes present in the BAM file, in string representation
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromBAM(string bamFile);

#endif