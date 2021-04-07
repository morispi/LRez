// gzindex is built above zran, as provided by Mark Adler.
// Full license follow.
/* zran.c -- example of zlib/gzip stream indexing and random access
 * Copyright (C) 2005, 2012 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
   Version 1.1  29 Sep 2012  Mark Adler */

/* Version History:
 1.0  29 May 2005  First version
 1.1  29 Sep 2012  Fix memory reallocation error
 */

/* Illustrate the use of Z_BLOCK, inflatePrime(), and inflateSetDictionary()
   for random access of a compressed file.  A file containing a zlib or gzip
   stream is provided on the command line.  The compressed stream is decoded in
   its entirety, and an index built with access points about every SPAN bytes
   in the uncompressed output.  The compressed file is left open, and can then
   be read randomly, having to decompress on the average SPAN/2 uncompressed
   bytes before getting to the desired block of data.

   An access point can be created at the start of any deflate block, by saving
   the starting file offset and bit of that block, and the 32K bytes of
   uncompressed data that precede that block.  Also the uncompressed offset of
   that block is saved to provide a referece for locating a desired starting
   point in the uncompressed stream.  build_index() works by decompressing the
   input zlib or gzip stream a block at a time, and at the end of each block
   deciding if enough uncompressed data has gone by to justify the creation of
   a new access point.  If so, that point is saved in a data structure that
   grows as needed to accommodate the points.

   To use the index, an offset in the uncompressed data is provided, for which
   the latest access point at or preceding that offset is located in the index.
   The input file is positioned to the specified location in the index, and if
   necessary the first few bits of the compressed data is read from the file.
   inflate is initialized with those bits and the 32K of uncompressed data, and
   the decompression then proceeds until the desired offset in the file is
   reached.  Then the decompression continues to read the desired uncompressed
   data from the file.

   Another approach would be to generate the index on demand.  In that case,
   requests for random access reads from the compressed data would try to use
   the index, but if a read far enough past the end of the index is required,
   then further index entries would be generated and added.

   There is some fair bit of overhead to starting inflation for the random
   access, mainly copying the 32K byte dictionary.  So if small pieces of the
   file are being accessed, it would make sense to implement a cache to hold
   some lookahead and avoid many calls to extract() for small lengths.

   Another way to build an index would be to use inflateCopy().  That would
   not be constrained to have access points at block boundaries, but requires
   more memory per access point, and also cannot be saved to file due to the
   use of pointers in the state.  The approach here allows for storage of the
   index in a file.
 */

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zlib.h"
#include <iostream>
#include <fstream>

#define local static

#define SPAN 1048576L       /* desired distance between access points */
#define WINSIZE 32768U      /* sliding window size */
#define CHUNK 16384         /* file input buffer size */

/* access point entry */
struct point {
    off_t out;          /* corresponding offset in uncompressed data */
    off_t in;           /* offset in input file of first full byte */
    int bits;           /* number of bits (1-7) from byte at in - 1, or 0 */
    unsigned char window[WINSIZE];  /* preceding 32K of uncompressed data */
};

/* access point list */
struct access {
    int have;           /* number of list entries filled in */
    int size;           /* number of list entries allocated */
    struct point *list; /* allocated list */
};

/** Deallocate an index built by build_index() 
    
    @param index index to deallocate
*/
void freeGzIndex(struct access *index);

/**
  Serialize the gzip index into a file.

  @param index index to serialize
  @param outputFile file where to store the index
*/
void serializeGzIndex(struct access* index, string outputFile);

/**
  Deserialize a gzip index from a file.

  @param index index to fill
  @param inputFile file where the index is stored
  @return the index populated with the contents of the file
*/
struct access* deserializeGzIndex(struct access* index, string inputFile);

/**
  Index a given gzip file.

  @param gzFile gzip file to build the index for
  @param span build access points every span bits of uncompressed output
  @param built structure where to store the index
  @return the number of access points on success (>=1), Z_MEM_ERROR for out of memory, 
  Z_DATA_ERROR for an error in the input file, or Z_ERRNO for a file read error.
  On success, *built points to the resulting index.
*/
int buildGzIndex(string, off_t span, struct access **built);

/**
  Index a given gzip file.

  @param FILE open on the desired gzip file to build the index for
  @param span build access points every span bits of uncompressed output
  @param built structure where to store the index
  @return the number of access points on success (>=1), Z_MEM_ERROR for out of memory, 
  Z_DATA_ERROR for an error in the input file, or Z_ERRNO for a file read error.
  On success, *built points to the resulting index.
*/
int buildGzIndex_Stream(FILE *in, off_t span, struct access **built);

/**
  Extract data starting from a given offset off the uncompressed file, using the index.

  @param gzFile gzip to extract data from
  @param index the index of the gzip file
  @param offset offset to start extracting data from
  @param buf buffer where to store the extracted data
  @param len number of bytes to extract
  @return bytes read or negative for error (Z_DATA_ERROR or Z_MEM_ERROR)
*/
int extract(string gzFile, struct access *index, off_t offset, unsigned char *buf, int len);

/**
  Extract data starting from a given offset off the uncompressed file, using the index.

  @param in FILE open on the desired gzip file
  @param index the index of the gzip file
  @param offset offset to start extracting data from
  @param buf buffer where to store the extracted data
  @param len number of bytes to extract
  @return bytes read or negative for error (Z_DATA_ERROR or Z_MEM_ERROR)
*/
int extract_Stream(FILE *in, struct access *index, off_t offset, unsigned char *buf, int len);

/**
  Extract the fastq format read at the specified offset.

  @param in FILE open on the desired gzip file
  @param index the index of the gzip file
  @param offset offset off the read to extract
  @return the extracted read in fastq format
*/
string extractFastqReadFromOffset(FILE* in, struct access* index, off_t offset);