#include "barcodesLoading.h"
#include "utils.h"

// Vectors used to translate Haplotagging barcodes to string
vector<string> Haplotagging_A;
vector<string> Haplotagging_B;
vector<string> Haplotagging_C;
vector<string> Haplotagging_D;

// Vector used to translate stLFR barcodes to string
vector<string> barcodes_stLFR;

void loadHaplotaggingBarcodes() {
    ifstream f;
    string line;
	f.open(path + "Haplotagging_A");
	if (!f.is_open()) {
		fprintf(stderr, "Unable to open barcodes list file barcodes/Haplotagging_A. Please provide an existing and valid file.\n");
		exit(EXIT_FAILURE);
	}

	while (getline(f, line)) {
		Haplotagging_A.push_back(line);
	}

	f.close();

	f.open(path + "Haplotagging_B");
	if (!f.is_open()) {
		fprintf(stderr, "Unable to open barcodes list file barcodes/Haplotagging_B. Please provide an existing and valid file.\n");
		exit(EXIT_FAILURE);
	}

	while (getline(f, line)) {
		Haplotagging_B.push_back(line);
	}

	f.close();

	f.open(path + "Haplotagging_C");
	if (!f.is_open()) {
		fprintf(stderr, "Unable to open barcodes list file barcodes/Haplotagging_C. Please provide an existing and valid file.\n");
		exit(EXIT_FAILURE);
	}

	while (getline(f, line)) {
		Haplotagging_C.push_back(line);
	}

	f.close();

	f.open(path + "Haplotagging_D");
	if (!f.is_open()) {
		fprintf(stderr, "Unable to open barcodes list file barcodes/Haplotagging_D. Please provide an existing and valid file.\n");
		exit(EXIT_FAILURE);
	}

	while (getline(f, line)) {
		Haplotagging_D.push_back(line);
	}

	f.close();

}

void loadstLFRBarcodes() {
    ifstream f;
    string line;

	f.open(path + "stLFR");
	if (!f.is_open()) {
		fprintf(stderr, "Unable to open barcodes list file barcodes/stLFR. Please provide an existing and valid file.\n");
		exit(EXIT_FAILURE);
	}

	while (getline(f, line)) {
		barcodes_stLFR.push_back(line);
	}

	f.close();
}