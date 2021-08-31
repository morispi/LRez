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

#include "gzIndex.h"

void serializeGzIndex(struct access* index, string outputFile) {
  ofstream out;
  out.open(outputFile, ios::out | ios::binary);
  if (!out.is_open()) {
    throw ios_base::failure("gzIndex: Unable to open file " + outputFile + ".");
  }

  unsigned char compData[WINSIZE];

  out.write((const char*) &index->have, sizeof(index->have));
  out.write((const char*) &index->size, sizeof(index->size));
  out.write((const char*) &index->maxOffset, sizeof(index->maxOffset));

  
  for (int i = 0; i < index->have; i++) {
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    deflateInit(&strm, Z_DEFAULT_COMPRESSION);
    strm.avail_in = WINSIZE;
    strm.next_in = index->list[i].window;
    strm.avail_out = WINSIZE;
    strm.next_out = compData;
    deflate(&strm, Z_FINISH);
    compData[WINSIZE - strm.avail_out] = '\0';

    out.write((const char*) &index->list[i].out, sizeof(index->list[i].out));
    out.write((const char*) &index->list[i].in, sizeof(index->list[i].in));
    out.write((const char*) &index->list[i].bits, sizeof(index->list[i].bits));
    auto size = WINSIZE - strm.avail_out;
    out.write((const char*) &size, sizeof(size));
    out.write((const char*) compData, WINSIZE - strm.avail_out);
    // out << endl;
    (void)deflateEnd(&strm);
    // out << compData << endl;
  }

  out.close();
}

struct access* deserializeGzIndex(struct access* index, string inputFile) {
  if (index == NULL) {
    index = (struct access*) malloc(sizeof(struct access));
  }

  ifstream in;
  in.open(inputFile, ios::in | ios::binary);
  if (!in.is_open()) {
    throw ios_base::failure("gzIndex: Untable to open gzip index for file " + inputFile + " for reading. Please make sure the gzip index file exists.");
  }

  in.read((char*) &index->have, sizeof(index->have));
  in.read((char*) &index->size, sizeof(index->size));
  in.read((char*) &index->maxOffset, sizeof(index->maxOffset));
  index->list = (struct point*) malloc(index->size * sizeof(struct point));

  for (int i = 0; i < index->have; i++) {
    in.read((char*) &index->list[i].out, sizeof(index->list[i].out));
    in.read((char*) &index->list[i].in, sizeof(index->list[i].in));
    in.read((char*) &index->list[i].bits, sizeof(index->list[i].bits));

    int size;
    in.read((char*) &size, sizeof(size));
    unsigned char buffer[WINSIZE];
    in.read((char*) buffer, size);
  
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    inflateInit(&strm);
    strm.avail_in = WINSIZE;
    strm.next_in = buffer;
    strm.avail_out = WINSIZE;
    unsigned char data[WINSIZE];
    strm.next_out = data;
    inflate(&strm, Z_FINISH);
    data[WINSIZE - strm.avail_out] = '\0';
    (void)inflateEnd(&strm);
    memcpy(index->list[i].window, data, WINSIZE - strm.avail_out);
    index->list[i].window[WINSIZE - strm.avail_out] = '\0';
  }


  return index;
}

void printAccess(struct access* index) {
  // cout << "have : " << index->have << endl;
  // cout << "size : " << index->size << endl;
    unsigned char tmp[WINSIZE];
    // tmp[0] = '\0';
    unsigned char tmp2[WINSIZE];
    // tmp2[0] = '\0';
  for (int i = 0; i < index->have; i++) {

      z_stream strm;
      strm.zalloc = Z_NULL;
      strm.zfree = Z_NULL;
      strm.opaque = Z_NULL;
      deflateInit(&strm, Z_DEFAULT_COMPRESSION);
      strm.avail_in = WINSIZE;
      strm.next_in = index->list[i].window;
      strm.avail_out = WINSIZE;
      strm.next_out = tmp;
      deflate(&strm, Z_FINISH);
      tmp[WINSIZE - strm.avail_out] = '\0';
      (void)deflateEnd(&strm);

      z_stream strm2;
      strm2.zalloc = Z_NULL;
      strm2.zfree = Z_NULL;
      strm2.opaque = Z_NULL;
      inflateInit(&strm2);
      strm2.avail_in = WINSIZE;
      strm2.next_in = tmp;
      strm2.avail_out = WINSIZE;
      strm2.next_out = tmp2;
      inflate(&strm2, Z_FINISH);
      tmp2[WINSIZE - strm2.avail_out] = '\0';
      (void)inflateEnd(&strm2);



    // cout << "out: " << index->list[i].out << endl;
    // cout << "in : " << index->list[i].in << endl;
    // cout << "bits : " << index->list[i].bits << endl;
    // cout << "window : " << tmp << endl;
    // cout << "window : " << tmp2 << endl;
  }
}

/* Add an entry to the access point list.  If out of memory, deallocate the
   existing list and return NULL. */
struct access *addpoint(struct access *index, int bits,
    off_t in, off_t out, unsigned left, unsigned char *window)
{
    struct point *next;

    /* if list is empty, create it (start with eight points) */
    if (index == NULL) {
        index = (struct access*) malloc(sizeof(struct access));
        if (index == NULL) return NULL;
        index->list = (struct point*) malloc(sizeof(struct point) << 3);
        if (index->list == NULL) {
            free(index);
            return NULL;
        }
        index->size = 8;
        index->have = 0;
    }

    /* if list is full, make it bigger */
    else if (index->have == index->size) {
        index->size <<= 1;
        next = (struct point*) realloc(index->list, sizeof(struct point) * index->size);
        if (next == NULL) {
            freeGzIndex(index);
            return NULL;
        }
        index->list = next;
    }

    /* fill in entry and increment how many we have */
    next = index->list + index->have;
    next->bits = bits;
    next->in = in;
    next->out = out;
    if (left)
        memcpy(next->window, window + WINSIZE - left, left);
    if (left < WINSIZE)
        memcpy(next->window + left, window, WINSIZE - left);
    index->have++;

    /* return list, possibly reallocated */
    return index;
}

/* Make one entire pass through the compressed stream and build an index, with
   access points about every span bytes of uncompressed output -- span is
   chosen to balance the speed of random access against the memory requirements
   of the list, about 32K bytes per access point.  Note that data after the end
   of the first zlib or gzip stream in the file is ignored.  build_index()
   returns the number of access points on success (>= 1), Z_MEM_ERROR for out
   of memory, Z_DATA_ERROR for an error in the input file, or Z_ERRNO for a
   file read error.  On success, *built points to the resulting index. */
int buildGzIndex_Stream(FILE *in, off_t span, struct access **built)
{
    int ret;
    off_t totin, totout;        /* our own total counters to avoid 4GB limit */
    off_t last;                 /* totout value of last access point */
    struct access *index;       /* access points being generated */
    z_stream strm;
    unsigned char input[CHUNK];
    unsigned char window[WINSIZE];

    /* initialize inflate */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2(&strm, 47);      /* automatic zlib or gzip decoding */
    if (ret != Z_OK)
        return ret;

    /* inflate the input, maintain a sliding window, and build an index -- this
       also validates the integrity of the compressed data using the check
       information at the end of the gzip or zlib stream */
    totin = totout = last = 0;
    index = NULL;               /* will be allocated by first addpoint() */
    strm.avail_out = 0;
    do {
        /* get some compressed data from input file */
        strm.avail_in = fread(input, 1, CHUNK, in);
        if (ferror(in)) {
            ret = Z_ERRNO;
            goto build_index_error;
        }
        if (strm.avail_in == 0) {
            ret = Z_DATA_ERROR;
            goto build_index_error;
        }
        strm.next_in = input;

        /* process all of that, or until end of stream */
        do {
            /* reset sliding window if necessary */
            if (strm.avail_out == 0) {
                strm.avail_out = WINSIZE;
                strm.next_out = window;
            }

            /* inflate until out of input, output, or at end of block --
               update the total input and output counters */
            totin += strm.avail_in;
            totout += strm.avail_out;
            ret = inflate(&strm, Z_BLOCK);      /* return at end of block */
            totin -= strm.avail_in;
            totout -= strm.avail_out;
            if (ret == Z_NEED_DICT)
                ret = Z_DATA_ERROR;
            if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR)
                goto build_index_error;
            if (ret == Z_STREAM_END)
                break;

            /* if at end of block, consider adding an index entry (note that if
               data_type indicates an end-of-block, then all of the
               uncompressed data from that block has been delivered, and none
               of the compressed data after that block has been consumed,
               except for up to seven bits) -- the totout == 0 provides an
               entry point after the zlib or gzip header, and assures that the
               index always has at least one access point; we avoid creating an
               access point after the last block by checking bit 6 of data_type
             */
            if ((strm.data_type & 128) && !(strm.data_type & 64) &&
                (totout == 0 || totout - last > span)) {
                index = addpoint(index, strm.data_type & 7, totin,
                                 totout, strm.avail_out, window);
                if (index == NULL) {
                    ret = Z_MEM_ERROR;
                    goto build_index_error;
                }
                last = totout;
            }
        } while (strm.avail_in != 0);
    } while (ret != Z_STREAM_END);



    /* clean up and return index (release unused entries in list) */
    (void)inflateEnd(&strm);
    index->list = (struct point*) realloc(index->list, sizeof(struct point) * index->have);
    index->size = index->have;
    // index->maxOffset = ftell(in);
    index->maxOffset = totout;
    *built = index;
    return index->size;

    /* return error */
  build_index_error:
    (void)inflateEnd(&strm);
    if (index != NULL)
        freeGzIndex(index);
    return ret;
}

int buildGzIndex(string gzFile, off_t span, struct access **built) {
  FILE *in;
  in = fopen(gzFile.c_str(), "rb");
  if (in == NULL) {
      throw ios_base::failure("gzIndex: could not open " + gzFile + " for reading. Please make sure the file exists.");

  }

  int res = buildGzIndex_Stream(in, span, built);

  fclose(in);
  return res;
}

/* Use the index to read len bytes from offset into buf, return bytes read or
   negative for error (Z_DATA_ERROR or Z_MEM_ERROR).  If data is requested past
   the end of the uncompressed data, then extract() will return a value less
   than len, indicating how much as actually read into buf.  This function
   should not return a data error unless the file was modified since the index
   was generated.  extract() may also return Z_ERRNO if there is an error on
   reading or seeking the input file. */
int extract_Stream(FILE* in, struct access *index, off_t offset,
                  unsigned char *buf, int len)
{
    int ret, skip;
    z_stream strm;
    struct point *here;
    unsigned char input[CHUNK];
    unsigned char discard[WINSIZE];

    /* proceed only if something reasonable to do */
    if (len < 0)
        return 0;

    /* find where in stream to start */
    here = index->list;
    ret = index->have;
    while (--ret && here[1].out <= offset)
        here++;

    /* initialize file and inflate state to start there */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2(&strm, -15);         /* raw inflate */
    if (ret != Z_OK)
        return ret;
    ret = fseeko(in, here->in - (here->bits ? 1 : 0), SEEK_SET);
    if (ret == -1)
        goto extract_ret;
    if (here->bits) {
        ret = getc(in);
        if (ret == -1) {
            ret = ferror(in) ? Z_ERRNO : Z_DATA_ERROR;
            goto extract_ret;
        }
        (void)inflatePrime(&strm, here->bits, ret >> (8 - here->bits));
    }
    (void)inflateSetDictionary(&strm, here->window, WINSIZE);

    /* skip uncompressed bytes until offset reached, then satisfy request */
    offset -= here->out;
    strm.avail_in = 0;
    skip = 1;                               /* while skipping to offset */
    do {
        /* define where to put uncompressed data, and how much */
        if (offset == 0 && skip) {          /* at offset now */
            strm.avail_out = len;
            strm.next_out = buf;
            skip = 0;                       /* only do this once */
        }
        if (offset > WINSIZE) {             /* skip WINSIZE bytes */
            strm.avail_out = WINSIZE;
            strm.next_out = discard;
            offset -= WINSIZE;
        }
        else if (offset != 0) {             /* last skip */
            strm.avail_out = (unsigned)offset;
            strm.next_out = discard;
            offset = 0;
        }

        /* uncompress until avail_out filled, or end of stream */
        do {
            if (strm.avail_in == 0) {
                strm.avail_in = fread(input, 1, CHUNK, in);
                if (ferror(in)) {
                    ret = Z_ERRNO;
                    goto extract_ret;
                }
                if (strm.avail_in == 0) {
                    ret = Z_DATA_ERROR;
                    goto extract_ret;
                }
                strm.next_in = input;
            }
            ret = inflate(&strm, Z_NO_FLUSH);       /* normal inflate */
            if (ret == Z_NEED_DICT)
                ret = Z_DATA_ERROR;
            if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR)
                goto extract_ret;
            if (ret == Z_STREAM_END)
                break;
        } while (strm.avail_out != 0);

        /* if reach end of stream, then don't keep trying to get more */
        if (ret == Z_STREAM_END)
            break;

        /* do until offset reached and requested data read, or stream ends */
    } while (skip);

    /* compute number of uncompressed bytes read after offset */
    ret = skip ? 0 : len - strm.avail_out;

    /* clean up and return bytes read or error */
  extract_ret:
    (void)inflateEnd(&strm);
    return ret;
}

int extract(string gzFile, struct access *index, off_t offset, unsigned char *buf, int len) {
  FILE *in;
  in = fopen(gzFile.c_str(), "rb");
  if (in == NULL) {
      throw ios_base::failure("gzIndex: could not open " + gzFile + " for reading. Please make sure the file exists.");
  }

  int res = extract_Stream(in, index, offset, buf, len);

  fclose(in);
  return res;
}

void freeGzIndex(struct access *index)
{
    if (index != NULL) {
        free(index->list);
        free(index);
    }
}

string extractFastqReadFromOffset(FILE* in, struct access* index, off_t offset) {
  unsigned char buf[CHUNK];
  extract_Stream(in, index, offset, buf, CHUNK);
  
  string str((const char*) buf);
  unsigned i = 0;
  unsigned nbNL = 0;
  while (nbNL < 4) {
    if (str[i] == '\n') {
      nbNL++;
    }
    i++;
  }

  return str.substr(0, i - 1);
}
