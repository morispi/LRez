# LRez

LRez provides a standalone tool allowing to work with barcoded linked-reads such as 10X Genomics data, as well as library allowing to easily use it in other projects.
LRez has different functionalities such as comparing regions pairs or contigs extremities to retrieve their common barcodes and extracting barcodes from given regions
of a BAM file, as well as indexing and querying both BAM and FASTQ files to quickly retrieve reads or alignments sharing a given barcode or list of barcodes.
In can thus be used in different applications, such as variant calling or scaffolding.

Requirements
--------------

  - A Linux based operating system.
  - g++, minimum version 5.5.0.
  - CMake, minimum version 2.8.2.
  - zlib, minimum version 1.2.11.
  - The Boost C++ library (https://www.boost.org/).
  
Installation
--------------

Clone the LRez repository, along with its submodules with:

  ```bash
  git clone --recursive https://github.com/morispi/LRez
  ```

Then run the install.sh script:

  ```bash
  ./install.sh
  ```

The installation script will build dependencies, the binary standalone in the `bin` folder, as well as the library, allowing to use LRez in other projects, in the `lib` folder.
  
Running LRez
--------------

### Usage

`LRez [SUBCOMMAND]`

where [SUBCOMMAND] can be one of the following:

  - compare:     Compute the number of common barcodes between pairs of regions, or between pairs of contigs' extremities
  - extract:     Extract the barcodes from a given region of a BAM file
  - index bam:   Index the offsets or occurrences positions of the barcodes contained in a BAM file
  - query bam:   Query the barcodes index to retrieve alignments in a BAM file, given a barcode or list of barcodes
  - index fastq: Index the offsets of the barcodes contained in a fastq file
  - query fastq: Query the barcodes index to retrieve alignments in a fastq file, given a barcode or list of barcodes

### Subcommands

A description of each subcommand as well as its options is given below.

#### Compare

LRez compare allows to compute the number of common barcodes between all possibles pairs of a given list of regions, or between a given contig's extremities and all other contigs' extremities.

      --bam STRING, -b STRING:      BAM file containing the alignments
      --index STRING, -i SRING:     Barcodes offsets index built with the index bam subcommand
      --region STRING, -r STRING:   File containing regions of interest in format chromosome:startPosition-endPosition
      --contig STRING, -c:          Contig of interest
      --size INT, -s INT:           Size of contigs' extremities to consider (optional, default: 1000) 
      --output STRING, -o STRING:   File where to output the results (optional, default: stdout)
      --userx, -u:                  Consider barcodes that only appear in the RX tag (optional, default: false)

#### Extract

LRez extract allows to extract the list of barcodes in a given region of a BAM file.

      --bam STRING, -b STRING:      BAM file to extract barcodes from
      --region STRING, -r STRING:   Region of interest in format chromosome:startPosition-endPosition
      --all, -a:                    Extract all barcodes
      --output STRING, -o STRING:   File where to output the extracted barcodes (optional, default: stdout)
      --userx, -u:                  Consider barcodes that only appear in the RX tag (optional, default: false)

#### Index BAM

LRez index bam allows to index the offsets or occurrences positions of the barcodes contained in a BAM file.

      --bam STRING, -b STRING:      BAM file to index
      --output STRING, -o STRING:   File where to store the index
      --offsets, -f:                Index the offsets of the barcodes in the BAM file
      --positions, -p:              Index the (chromosome, begPosition) occurrences positions of the barcodes
      --primary, -r:                Only index barcodes that appear in a primary alignment (optional, default: false)
      --quality INT, -q INT:        Only index barcodes that appear in an alignment of quality higher than this number (optional, default: 0)
      --userx, -u:                  Consider barcodes that only appear in the RX tag (optional, default: false)

#### Query BAM

LRez query bam allows to query a barcodes index and a BAM file to retrieve alignments containing the query barcodes.

      --bam STRING, -b STRING:      BAM file to search
      --index STRING, -i STRING:    Barcodes offsets index, built with the index bam subcommand.
      --list STRING, -l STRING:     File containing a list of barcodes to search in the BAM / index
      --output STRING, -o STRING:   File where to output the extracted alignments (optional, default: stdout)
      --userx, -u:                  Consider barcodes that only appear in the RX tag (optional, default: false)

#### Index fastq

LRez index fastq allows to index the offsets of the barcodes contained in a fastq file.

      --fastq STRING, -f STRING:    Fastq file to index
      --output STRING, -o STRING:   File where to store the index
      --gzip, -g:                   Fastq file is gzipped (optional, default: false)
      --userx, -u:                  Consider barcodes that only appear in the RX tag (optional, default: false)

#### Query fastq

LRez query fastq allows to query a barcodes index and a fastq file to retrieve alignments containing the query barcodes.

      --fastq STRING, -f STRING:    Fastq file to search
      --index STRING, -i STRING:    Barcodes index, built with the index fastq subcommand
      --query STRING, -q STRING:    Query barcode to search in the fastq file and the index
      --list STRING, -l STRING:     File containing a list of barcodes to search in the fastq file and the index
      --output STRING, -o STRING:   File where to output the extracted reads (optional, default: stdout)
      --gzip, -g:                   Fastq file is gzipped (optional, default: false)
      --userx, -u:                  Consider barcodes that only appear in the RX tag (optional, default: false)


Documentation
--------------

Complete documentation of the different functions is provided in the `.h` files of the `include/` and `include/subcommands/` folders. Proper documentation will be provided on this page later.

Notes
--------------

LRez has been developed and tested on x86-64 GNU/Linux.          
Support for any other platform has not been tested.

Authors
--------------

Pierre Morisse, Claire Lemaitre and Fabrice Legeai.

Reference
--------------

Pierre Morisse, Claire Lemaitre, Fabrice Legeai. LRez: C++ API and toolkit for analyzing and managing Linked-Reads data. https://arxiv.org/abs/2103.14419

Contact
--------------

You can report problems and bugs to pierre[dot]morisse[at]inria[dot]fr
