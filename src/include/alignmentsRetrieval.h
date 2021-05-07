#ifndef __LREZ_ALIGNMENTS_RETRIEVAL__
#define __LREZ_ALIGNMENTS_RETRIEVAL__

#include "utils.h"
#include "indexManagementBam.h"

/**
	Retrieve the alignment starting at the given position and containing the given barcode.

	@param alignment pointer to the variable where to store the alignment if it is found
	@param reader BamReader open on the desired BAM file
	@param position position at which the alignment must start
	@param b barcode of interest, in binary representation
	@return the alignment of interest
*/
BamAlignment retrieveAlignmentWithBarcode(BamAlignment* alignment, BamReader& reader, int64_t position, barcode b);


/**
	Retrieve all the alignments that share a given barcode, in a given BAM file.
	A barcode index of the BAM must be available.

	@param reader BamReader open on the desired BAM file
	@param BarcodesOffsetsIndex barcode offsets index of the BAM file
	@param b barcode of interest, in binary representation
	@return a vector containing all the alignments sharing the barcode of interest
*/
vector<BamAlignment> retrieveAlignmentsWithBarcodeBits_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, barcode b);


/**
	Retrieve all the alignments that share a given barcode, in a given BAM file.
	A barcode index of the BAM must be available.

	@param reader BamReader open on the desired BAM file
	@param BarcodesOffsetsIndex barcode offsets index of the BAM file
	@param bc barcode of interest
	@return a vector containing all the alignments sharing the barcode of interest
*/
vector<BamAlignment> retrieveAlignmentsWithBarcode_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string bc);

/**
	Retrieve all the alignments that share a given barcode, in a given BAM file.
	A barcode index of the BAM must be available.

	@param bamFile BAM file to extract alignments from
	@param BarcodesOffsetsIndex barcode offsets index of the BAM file
	@param bc barcode of interest
	@return a vector containing all the alignments sharing the barcode of interest
*/
vector<BamAlignment> retrieveAlignmentsWithBarcode(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string bc);

/**
	Retrieve all the alignments that share a given barcode, in a given BAM file.
	A barcode index of the BAM must be available.

	@param bamFile BAM file to extract alignments from
	@param BarcodesOffsetsIndex barcode offsets index of the BAM file
	@param bc barcode of interest, in binary representation
	@return a vector containing all the alignments sharing the barcode of interest
*/
vector<BamAlignment> retrieveAlignmentsWithBarcodeBits(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, barcode bc);

/**
	Retrieve all the alignments that share a barcode appearing in the list of a given file

	@param in stream open on the BAM fastq file
	@param BarcodesOffsetsIndex barcode offsets index of the BAM file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line 
	@return a vector containing all the alignments with the barcodes of interest
*/
vector<BamAlignment> retrieveAlignmentsWithBarcodes_BamReader(BamReader& reader, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string barcodesList);

/**
	Retrieve all the alignments that share a barcode appearing in the list of a given file

	@param bamFile BAM file to extract alignments from
	@param BarcodesOffsetsIndex barcode offsets index of the BAM file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line 
	@return a vector containing all the alignments with the barcodes of interest
*/
vector<BamAlignment> retrieveAlignmentsWithBarcodes(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string barcodesList);

#endif