#ifndef __LREZ_UTILS__
#define __LREZ_UTILS__

#include <string>
#include <vector>
#include "api/BamReader.h"
#include "api/BamIndex.h"
#include "api/BamAux.h"
#include "robin_hood.h"
#include <sstream>
#include <regex>
#include "barcodesLoading.h"

#define BXTAG "BX:Z"

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

using namespace std;
using namespace BamTools;

typedef vector<bool> barcode;

/**
 Supported sequencing technologies
*/
enum SequencingTechnology {Undefined = 0, TenX, Haplotagging, TELLSeq, stLFR};

extern SequencingTechnology techno;

// Vectors used to translate Haplotagging barcodes to string
extern vector<string> Haplotagging_A;
extern vector<string> Haplotagging_B;
extern vector<string> Haplotagging_C;
extern vector<string> Haplotagging_D;

// Vector used to translate stLFR barcodes to string
extern vector<string> barcodes_stLFR;

/**
	Determine the sequencing technology the barcode originates from.
	This function is compatible with 10x Genomics, Haplotagging, TELL-SEq and stLFR.
	Barcodes that do not come from these technologies and are not represented as a suite of nucleotides will cause the program to exit.

	@param barcode the barcode to determine sequencing technology from
	@return The SequencingTechnology enum field corresponding to the sequencing technology
*/
SequencingTechnology determineSequencingTechnology(const string& barcode);

/**
	Retrieve the nucleotides content of the barcode.
	This function is used to translate barcodes represented as a suite of integers (as in stLFT and Haplotagging) into nucleotides barcodes.
	The function takes care of determining the employed sequencing technoly.

	@param barcode the barcode to retrieve nucleotides for
	@return the barcode in nucleotides representation
*/
string retrieveNucleotidesContent(const string& barcode);

/**
	Check whether a barcode is valid or not.
	A barcode is considered as valid if it is not empty, if it does not contain any "N" for 10x and TELL-Seq,
	if it is not "0_0_0" for stLFR data, and does not contain a "00" substring for Haplotagging data.
	The function takes care of determining the employed sequencing technoly.

	@param barcode the barcode to verify
	@return true if the barcode is valid, false otherwise
*/
bool isValidBarcode(const string& barcode);

/**
	Translate a string to a barcode in 2 bits per nucleotide format.
	The function takes care of determining the employed sequencing technoly, and of retrieving the nucleotides contents of barcodes represented as a suite of integers.

	@param str string to convert
	@return the barcode in binary representation
*/
barcode stringToBarcode(const string& str);

/**
	Split a string according to a delimiter.

	@param s string to split
	@param delimiter delimiter 
	@return a vector containing the splits of the string
*/
vector<string> splitString(string s, string delimiter);

/**
	Translate a string to a BamRegion.

	@param reader BamReader open on the desired BAM file
	@param s string to translate, formatted as chromosome:startPosition-endPosition
	@return the BamRegion coressponding to the string s
*/
BamRegion stringToBamRegion(BamReader& reader, string s);

/**
	Translate a BamAlignment to a SAM-like string.

	@param a BamAlignment to translate
	@param m_references vector containing the information (name and length) about reference sequences
	@return a SAM-like string summarizing the information of a
*/
string convertToSam(const BamAlignment& a, RefVector m_references);

#endif
