#ifndef __LREZ_BARCODES_EXTRACTION__
#define __LREZ_BARCODES_EXTRACTION__

#include "utils.h"
#include "robin_hood.h"

/**
	Extract the barcodes of a given region, in binary representation.

	@param reader BamReader open on the desired BAM file
	@param region region of interest
	@return a set containing the barcodes present in the region
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion_BamReader(BamReader& reader, string region);

/**
	Extract the barcodes of a given region, in string representation.

	@param reader BamReader open on the desired BAM file
	@param region region of interest
	@return a set containing the barcodes present in the region
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion_BamReader(BamReader& reader, string region);

/**
	Extract the barcodes of a given region, in binary representation.

	@param reader bamFile BAM file to extract barcodes from
	@param region region of interest
	@return a set containing the barcodes present in the region, in binary representation
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromRegion(string bamFile, string region);


/**
	Extract the barcodes of a given region, in string representation.

	@param reader bamFile BAM file to extract barcodes from
	@param region region of interest
	@return a set containing the barcodes present in the region, in binary string
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromRegion(string bamFile, string region);

/**
	Extract all the barcodes of a given BAM file, in binary representation.

	@param bamFile bamFile BAM file to extract barcodes from
	@return a set containing the barcodes present in the BAM file, in binary representation
*/
robin_hood::unordered_set<barcode> extractBarcodesBitsFromBAM(string bamFile);

/**
	Extract all the barcodes of a given BAM file, in string representation.

	@param bamFile bamFile BAM file to extract barcodes from
	@return a set containing the barcodes present in the BAM file, in string representation
*/
robin_hood::unordered_set<string> extractBarcodesSeqsFromBAM(string bamFile);

/**
	Extract the barcodes of a given region, in binary representation, including duplicates.

	@param reader BamReader open on the desired BAM file
	@param region region of interest
	@return a vector containing the barcodes present in the region
*/
vector<barcode> extractBarcodesBitsFromRegionWithDuplicates_BamReader(BamReader& reader, string region);

/**
	Extract the barcodes of a given region, in string representation, including duplicates.

	@param reader BamReader open on the desired BAM file
	@param region region of interest
	@return a vector containing the barcodes present in the region
*/
vector<string> extractBarcodesSeqsFromRegionWithDuplicates_BamReader(BamReader& reader, string region);

/**
	Extract the barcodes of a given region, in binary representation, including duplicates.

	@param reader bamFile BAM file to extract barcodes from
	@param region region of interest
	@return a vector containing the barcodes present in the region, in binary representation
*/
vector<barcode> extractBarcodesBitsFromRegionWithDuplicates(string bamFile, string region);


/**
	Extract the barcodes of a given region, in string representation, including duplicates.

	@param reader bamFile BAM file to extract barcodes from
	@param region region of interest
	@return a vector containing the barcodes present in the region, in binary string
*/
vector<string> extractBarcodesSeqsFromRegionWithDuplicates(string bamFile, string region);

/**
	Extract all the barcodes of a given BAM file, in binary representation, including duplicates.

	@param bamFile bamFile BAM file to extract barcodes from
	@return a vector containing the barcodes present in the BAM file, in binary representation
*/
vector<barcode> extractBarcodesBitsFromBAMWithDuplicates(string bamFile);

/**
	Extract all the barcodes of a given BAM file, in string representation, including duplicates.

	@param bamFile bamFile BAM file to extract barcodes from
	@return a vector containing the barcodes present in the BAM file, in string representation
*/
vector<string> extractBarcodesSeqsFromBAMWithDuplicates(string bamFile);

#endif