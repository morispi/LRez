#ifndef __LREZ_READS_RETRIEVAL__
#define __LREZ_READS_RETRIEVAL__

#include "utils.h"
#include "indexManagementFastq.h"


/**
	Retrieve all the reads that share a given barcode, in a given fastq file.

	@param in stream open on the desired fastq file
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest, in binary representation
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcodeBits_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, barcode bc);

/**
	Retrieve all the reads that share a given barcode, in a given fastq file.

	@param in stream open on the desired fastq file
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcode_Stream(ifstream& reader, BarcodesIndex& BarcodesIndex, string bc);

/**
	Retrieve all the reads that share a given barcode, in a given fastq file.

	@param fastqFile fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest, in binary representation
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcodeBits(string fastqFile, BarcodesIndex& BarcodesIndex, barcode bc);

/**
	Retrieve all the reads that share a given barcode, in a given fastq file.

	@param fastqFile fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcode(string fastqFile, BarcodesIndex& BarcodesIndex, string bc);

/**
	Retrieve all the reads that share a barcode appearing in the list of a given file

	@param in stream open on the desired fastq file
	@param BarcodesIndex barcode index of the fastq file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line 
	@return a vector containing all the reads with the barcodes of interest
*/
vector<string> retrieveReadsWithBarcodes_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, string barcodesList);

/**
	Retrieve all the reads that share a barcode appearing in the list of a given file

	@param fastqFile fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line 
	@return a vector containing all the reads with the barcodes of interest
*/
vector<string> retrieveReadsWithBarcodes(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList);

#endif