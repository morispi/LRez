import pysam
import sys
import os

bamFile = pysam.AlignmentFile(sys.argv[1], "rb");
outFile = pysam.AlignmentFile(os.path.splitext(sys.argv[1])[0] + "_barcodes_extracted" + os.path.splitext(sys.argv[1])[1], "wb", template=bamFile)

iter = bamFile.fetch()
for al in iter:
        t = al.query_name.split("#")
        al.query_name = t[0]
        al.set_tag("BX", t[1])
        outFile.write(al)

bamFile.close()
outFile.close()
