#ifndef __LREZ_READS_RETRIEVAL__
#define __LREZ_READS_RETRIEVAL__

#include "utils.h"
#include "indexManagementFastq.h"
#include "gzIndex.h"

/**
	Retrieve all the reads that share a given barcode, in a given fastq file.

	@param in stream open on the desired fastq file
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest, in binary representation
	@throws runtime_error thrown if a given offset of the file where the barcode appears could not be seeked
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcodeBits_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, barcode bc);

/**
	Retrieve all the reads that share a given barcode, in a given fastq file.

	@param fastqFile fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest, in binary representation
	@throws ios_base::failure thrown if fastqFile could not be open
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcodeBits(string fastqFile, BarcodesIndex& BarcodesIndex, barcode bc);

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
	@param bc barcode of interest
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcode(string fastqFile, BarcodesIndex& BarcodesIndex, string bc);

/**
	Retrieve all the reads of a fastq file that share a barcode appearing in the list of a given file

	@param in stream open on the desired fastq file
	@param BarcodesIndex barcode index of the fastq file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line
	@throws ios_base::failure thrown if barcodesList could not be open
	@return a vector containing all the reads with the barcodes of interest
*/
vector<string> retrieveReadsWithBarcodes_Stream(ifstream& in, BarcodesIndex& BarcodesIndex, string barcodesList);

/**
	Retrieve all the reads of a fastq file that share a barcode appearing in the list of a given file

	@param fastqFile fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line
	@param nbThreads number of threads to use, set to 1 by default
	@throws ios_base::failure thrown if fastqFile or barcodesList could not be open
	@return a vector containing all the reads with the barcodes of interest
*/
vector<string> retrieveReadsWithBarcodes(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList, unsigned nbThreads = 1);


/**
	Retrieve all the reads that share a given barcode, in a given gzipped fastq file, using its index.

	@param in FILE open on the desired gzipped fastq file
	@param gzIndex index of the gzipped fastq file
	@param BarcodesIndex barcode index of the gzipped fastq file
	@param bc barcode of interest, in binary representation
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcodeBits_Gzip_Stream_Index(FILE* in, struct access* gzIndex, BarcodesIndex& barcodexIndex, barcode bc);

/**
	Retrieve all the reads that share a given barcode, in a given gzipped fastq file, using its index.

	@param fastqFile gzipped fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest, in binary representation
	@throws ios_base::failure thrown if fastqFile could not be open
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcodeBits_Gzip(string fastqFile, BarcodesIndex& barcodexIndex, barcode bc);


/**
	Retrieve all the reads that share a given barcode, in a given gzipped fastq file, using its index.

	@param in FILE open on the desired gzipped fastq file
	@param gzIndex index of the gzipped fastq file
	@param BarcodesIndex barcode index of the gzipped fastq file
	@param bc barcode of interest
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcode_Gzip_Stream_Index(FILE* in, struct access* gzIndex, BarcodesIndex& BarcodesIndex, string barcode);


/**
	Retrieve all the reads that share a given barcode, in a given gzipped fastq file, using its index.

	@param fastqFile gzipped fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param bc barcode of interest
	@return a vector containing all the reads with the barcode of interest
*/
vector<string> retrieveReadsWithBarcode_Gzip(string fastqFile, BarcodesIndex& BarcodesIndex, string barcode);

/**
	Retrieve all the reads of a gzipped fastq file that share a barcode appearing in the list of a given file

	@param in FILE open on the desired gzipped fastq file
	@param gzIndex index of the gzipped fastq file
	@param BarcodesIndex barcode index of the fastq file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line 
	@throws ios_base::failure thrown if barcodesList could not be open
	@return a vector containing all the reads with the barcodes of interest
*/
vector<string> retrieveReadsWithBarcodes_Gzip_Stream_Index(FILE* in, struct access* gzIndex, BarcodesIndex& BarcodesIndex, string barcodesList);

/**
	Retrieve all the reads of a fastq file that share a barcode appearing in the list of a given file

	@param fastqFile fastq file to extract reads from
	@param BarcodesIndex barcode index of the fastq file
	@param barcodesList file containing the list of barcodes of interest, with one barcode per line 
	@param nbThreads number of threads to use, set to 1 by default
	@throws ios_base::failure thrown if fastqFile or barcodesList could not be open
	@return a vector containing all the reads with the barcodes of interest
*/
vector<string> retrieveReadsWithBarcodes_Gzip(string fastqFile, BarcodesIndex& BarcodesIndex, string barcodesList, unsigned nbThreads = 1);



#endif