#ifndef __LREZ_UTILS__
#define __LREZ_UTILS__

#include <string>
#include <vector>
#include "api/BamReader.h"
#include "api/BamIndex.h"
#include "api/BamAux.h"
#include "robin_hood.h"
#include <sstream>

//TODO: See how to define barcodes
#define BARCODE_SIZE 16
#define BXTAG "BX:Z"
// BX tag is not always there
#define RXTAG "RX:Z"

#define no_argument 0
#define required_argument 1 
#define optional_argument 2

using namespace std;
using namespace BamTools;

typedef uint32_t barcode;

extern bool CONSIDER_RX;


/**
	Translate a string to a barcode.

	@param str string to convert
*/
barcode stringToBarcode(const string& str);

/**
	Translate a barcode to a string.

	@param b barcode to convert
*/
string barcodeToString(barcode b);

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
