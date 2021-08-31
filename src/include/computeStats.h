#ifndef __LREZ_COMPUTE_STATS__
#define __LREZ_COMPUTE_STATS__

#include "utils.h"
#include <map>

/**
	Structure recording global stats.
	nbBarcodes: Overall number of barcodes
	nbMappedReads: Overrall number of mapped reads
	readsPerBarcode: List of the counts of number of reads per barcode
*/
struct Stats {
	unsigned nbBarcodes;
	unsigned nbMappedReads;
	vector<unsigned> readsPerBarcode;
};


/**
	Extract the number of barcodes from a random list of regions.
	@param reader BamReader open on the desired BAM file
	@param regionsList the list of all regions of the reference genome
	@param numberOfRegions number of regions to extract barcodes from
	@return A sorted list containing the number of barcodes of the randomly selected regions
*/
vector<unsigned> extractBarcodesCountsPerRegion(BamReader& reader, vector<string> regionsList, unsigned numberOfRegions);

/**
	Extract the number of common barcodes between random adjacent regions.
	@param reader BamReader open on the desired BAM file
	@param regionsList the list of all regions of the reference genome
	@param numberOfRegions number of adjacent windows to consider
	@return A sorted list containing the number of common barcodes between the randomly selected adjacent regions
*/
vector<unsigned> extractCommonBarcodesCounts(BamReader& reader, vector<string> regionsList, unsigned numberOfRegions);

/**
	Extract the number of barcodes from a random list of regions and the number of common barcodes between adjacent regions of these regions.
	@param bamFile BAM file to extract stats from
	@param regionsList the list of all regions of the reference genome
	@param numberOfRegions number of adjacent windows to consider
	@param nbThreads number of threads to use, set to 1 by default
	@throws ios_base::failure thrown if bamFile or its associated index could not be open
	@return A pair composed of sorted list containing the number of barcodes of the randomly selected regions and a sorted list containing the number of common barcodes between adjacent windows of the selected windows
*/
pair<vector<unsigned>, vector<unsigned>> extractBarcodesAndCommonBarcodesCounts(string bamFile, vector<string> regionsList, unsigned numberOfRegions, unsigned nbThreads = 1);

/**
	Extract global stats (number of barcodes, number of mapped reads, number of reads for each barcode).
	@param bamFile BAM file to extract stats from
	@param nbThreads number of threads to use, set to 1 by default
	@throws ios_base::failure thrown if bamFile or its associated index, or if regions could not be open
	@throws runtime_error thrown if bamFile could not be seeked
	@return A structure describing the global stats of the BAM file
*/
Stats extractGlobalStats(string bamFile, unsigned nbThreads = 1);

#endif
