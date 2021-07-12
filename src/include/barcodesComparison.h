#ifndef __LREZ_BARCODES_COMPARISON__
#define __LREZ_BARCODES_COMPARISON__

#include "utils.h"
#include "indexManagementBam.h"

/**
	Structure defining a hash function, to hash a pair of any kind.
	Required to use a pair<string, string> as key in the maps used for the compareRegions and compareContig functions.
*/
struct hashPairs { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

// struct hashPairsPairs { 
//     template <class T1, class T2> 
//     size_t operator()(const pair<T1, T2>& p) const
//     { 
//         auto hash1 = p.first.first ^ p.first.second; 
//         auto hash2 = p.second.first ^ p.second.second; 
//         return hash1 ^ hash2; 
//     } 
// }; 

/**
	Count the number of common barcodes between two sets of barcodes.

	@param barcodes1 first set of barcodes
	@param barcodes2 second set of barcodes
	@return the number of common barcodes between the two sets
*/
unsigned countCommonBarcodes(robin_hood::unordered_set<barcode> barcodes1, robin_hood::unordered_set<barcode> barcodes2);

/**
	Compute the number of barcodes that are shared between all possible pair of regions contained in a map.

	@param regionsBarcodes map associating every region of interest to its set of barcodes
	@return a map associating all possible pairs of regions to their number of common barcodes 
*/
robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> computePairwiseCommonBarcounts(robin_hood::unordered_map<string, robin_hood::unordered_set<barcode>> regionsBarcodes);

/**
	Compute the number of barcodes that are shared between all possible pair of regions contained in a file.

	@param reader bamFile BAM file to compute common barcodes from
	@param regions file containing a list of regions formatted as chromosome:startPosition-endPosition
	@return a map associating all possible pairs of regions to their number of common barcodes 
*/
robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareRegions(string bamFile, string regions);

/**
	Compute the number of barcodes that are shared between a given contig's extremity and other contigs' extremities.

	@param counts map associating pairs of contigs' extremities to their number of common barcodes
	@param BarcodesOffsetsIndex barcode index associating barcodes to the set of BamRegion they appear in
	@param reader BamReader open on the desired BAM file
	@param rv vector containing the information (name and length) about reference sequences
	@param size size of the contigs extremities to consider
	@param qRegion region of the contig extremity of interest 
*/
void computeCommonBarcodesCounts(robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs>& counts, BarcodesOffsetsIndex& BarcodesOffsetsIndex, BamReader& reader, const RefVector& rv, int size, string& qRegion);

/**
	Compare the barcodes of a given contig's extremities to the barcodes of other contigs' extremities.

	@param reader BamReader open on the desired BAM file
	@param BarcodesOffsetsIndex barcode index associating barcodes to the set of BamRegion they appear in
	@param contig name of the contig of interest
	@param size size of the contigs extremities to consider
	@return a map associating pairs of contigs' extremities to their number of common barcodes
*/
robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContig_BamReader(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contig, int size);

/**
	Compare the barcodes of a given contig's extremities to the barcodes of other contigs' extremities.

	@param bamFile BAM file to compare contigs' extremities barcodes from
	@param BarcodesOffsetsIndex barcode index associating barcodes to the set of BamRegion they appear in
	@param contig name of the contig of interest
	@param size size of the contigs extremities to consider
	@return a map associating pairs of contigs' extremities to their number of common barcodes
*/
robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContig(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contig, int size);

/**
	Compare the barcodes of a given list of contigs' extremities to the barcodes of other contigs' extremities.

	@param bamFile BAM file to compare contigs' extremities barcodes from
	@param BarcodesOffsetsIndex barcode index associating barcodes to the set of BamRegion they appear in
	@param contigs file containing a list of contigs of interest
	@param size size of the contigs extremities to consider
	@return a map associating pairs of contigs' extremities to their number of common barcodes
*/
robin_hood::unordered_map<pair<string, string>, unsigned, hashPairs> compareContigs(string bamFile, BarcodesOffsetsIndex& BarcodesOffsetsIndex, string contigs, int size);

#endif