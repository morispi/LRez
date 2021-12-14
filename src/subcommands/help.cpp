#include "subcommands/help.h"

void subcommandHelp(std::string subcommand) {
	if (subcommand == "global") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez allows to work with barcoded Linked-Reads, and offers various barcode management functionalities.\n");
		printf("\n");

		printf("USAGE:\n");
		printf("\tLRez [SUBCOMMAND]\n");
		printf("\n");

		printf("SUBCOMMANDS:\n");
		printf("\tcompare\t\t Compute the number of common barcodes between pairs of regions, or between pairs of contigs' extremities\n");
		printf("\textract\t\t Extract the barcodes from a given region of a BAM file\n");
		printf("\tstats\t\t Retrieve general stats from a BAM file\n");
		printf("\tindex bam\t Index the offsets or occurrences positions of the barcodes contained in a BAM file\n");
		printf("\tquery bam\t Query the barcodes index to retrieve alignments in a BAM file, given a barcode or list of barcodes\n");
		printf("\tindex fastq\t Index the offsets of the barcodes contained in a fastq file\n");
		printf("\tquery fastq\t Query the barcodes index to retrieve alignments in a fastq file, given a barcode or list of barcodes\n");
	} else if (subcommand == "compare") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez compare allows to compute the number of common barcodes between all possibles pairs of a given list of regions, or between a given contig's extremities and all other contigs' extremities\n");	
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez compare [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-b, --bam\t BAM file containing the alignments\n");
		printf("\t-i, --index\t Barcodes offsets index built with the index bam subcommand\n");
		printf("\t-r, --regions\t File containing regions of interest in format chromosome:startPosition-endPosition\n");
		printf("\t-c, --contig\t Contig of interest\n");
		printf("\t-C, --contigs\t File containing a list of contigs of interest\n");
		printf("\t-s, --size\t Size of contigs' extremities to consider (optional, default: 1000)\n");
		printf("\t-o, --output\t File where to output the results (optional, default: stdout)\n");
		printf("\t-t, --threads\t Number of threads to use when comparing a list of contigs (optional, default: 1)\n");
	} else if (subcommand == "extract") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez extract allows to extract the list of barcodes in a given region of a BAM file.\n");	
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez extract [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-b, --bam\t\t BAM file to extract barcodes from\n");
		printf("\t-r, --region\t\t Region of interest in format chromosome:startPosition-endPosition\n");
		printf("\t-a, --all\t\t Extract all barcodes\n");
		printf("\t-o, --output\t\t File where to output the extracted barcodes (optional, default: stdout)\n");
		printf("\t-d, --duplicates\t Include duplicate barcodes (optional, default: false)\n");
		printf("\t-t, --threads\t Number of threads to use when comparing a list of contigs (optional, default: 1)\n");
	} else if (subcommand == "index bam") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez index bam allows to index the offsets or occurrences positions of the barcodes contained in a BAM file.\n");
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez index bam [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-b, --bam\t BAM file to index\n");
		printf("\t-o, --output\t File where to store the index\n");
		printf("\t-f, --offsets\t Index the offsets of the barcodes in the BAM file\n");
		printf("\t-p, --positions\t Index the (chromosome, begPosition) occurrences positions of the barcodes\n");
		printf("\t-r, --primary\t Only index barcodes that appear in a primary alignment (optional, default: false)\n");
		printf("\t-q, --quality\t Only index barcodes that appear in an alignment of quality higher than this number (optional, default: 0)\n");
		printf("\t-t, --threads\t Number of threads to use to build the index (optional, default: 1)\n");
	} else if (subcommand == "query bam") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez query bam allows to query a barcodes index and a BAM file to retrieve alignments containing the query barcodes.\n");
		printf("Matched alignments are returned in SAM format.\n");	
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez query bam [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-b, --bam\t BAM file to search\n");
		printf("\t-i, --index\t Barcodes offsets index, built with the index bam subcommand.\n");
		printf("\t-q, --query\t Query barcode to search in the BAM / index\n");
		printf("\t-l, --list\t File containing a list of barcodes to search in the BAM / index\n");
		printf("\t-o, --output\t File where to output the extracted alignments (optional, default: stdout)\n");
		printf("\t-H, --header\t Output SAM header (optional, default: false)\n");
		printf("\t-t, --threads\t Number of threads to use when querying with a list of barcodes (optional, default: 1)\n");
	} else if (subcommand == "index fastq") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez index fastq allows to index the offsets of the barcodes contained in a fastq file.\n");
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez index fastq [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-f, --fastq\t Fastq file to index\n");
		printf("\t-o, --output\t File where to store the index\n");
		printf("\t-g, --gzip\t Fastq file is gzipped (optional, default: false)\n");
		printf("\t-t, --threads\t Number of threads to use to build the index (optional, default: 1)\n");
	} else if (subcommand == "query fastq") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez query fastq allows to query a barcodes index and a fastq file to retrieve alignments containing the query barcodes.\n");
		printf("Matched alignments are returned in SAM format.\n");	
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez query fastq [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-f, --fastq\t Fastq file to search\n");
		printf("\t-i, --index\t Barcodes index, built with the index fastq subcommand\n");
		printf("\t-q, --query\t Query barcode to search in the fastq file and the index\n");
		printf("\t-l, --list\t File containing a list of barcodes to search in the fastq file and the index\n");
		printf("\t-c, --collectionOfLists\t File of files (FOF) e.g. file containing files' names of lists of barcodes to search in the fastq file and the index\n");
		printf("\t-o, --output\t File where to output the extracted reads (optional, default: stdout)\n");
		printf("\t-g, --gzip\t Fastq file is gzipped (optional, default: false)\n");
		printf("\t-t, --threads\t Number of threads to use when querying with a list of barcodes (optional, default: 1)\n");
	} else if (subcommand == "stats") {
		printf("%s\n", VERSION);
		printf("Pierre Morisse <pierre.morisse@inria.fr>\n");
		printf("LRez stats allows to retrieve general stats from a BAM file.\n");
		printf("\n");
		
		printf("USAGE:\n");
		printf("\tLRez stats [ARGS]\n");
		printf("\n");

		printf("ARGS:\n");
		printf("\t-b, --bam\t BAM file to retrieve stats from\n");
		printf("\t-r, --regions\t Number of regions to consider to define stats (optional, default: 1000)\n");
		printf("\t-s, --size\t Size of the regions to consider (optional, default: 1000) \n");
		printf("\t-o, --output\t File where to output the extracted reads (optional, default: stdout)\n");
		printf("\t-t, --threads\t Number of threads to use when comparing a list of contigs (optional, default: 1)\n");
	}

	exit(EXIT_FAILURE);
}
