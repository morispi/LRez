#ifndef __LREZ_BARCODES_LOADING__
#define __LREZ_BARCODES_LOADING__

#include "utils.h"

using namespace std;

// Path were the barcodes files are stored
extern string path;

/**
	Populate the Haplotagging_[ABCD] vectors to Haplotagging barcodes can be coverted to nucleotides.
*/
void loadHaplotaggingBarcodes();

/**
	Populate the barcodes_stLFR vector to Haplotagging barcodes can be coverted to nucleotides.
*/
void loadstLFRBarcodes();

#endif