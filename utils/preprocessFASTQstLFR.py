import sys
import re
import csv
import os

r = open(sys.argv[1]);
o = open(os.path.splitext(sys.argv[1])[0] + "_barcodes_in_headers" + os.path.splitext(sys.argv[1])[1], "w");

header = r.readline()[:-1];
while header != '':
	line = r.readline();
	data = line;
	line = r.readline();
	data += line;
	line = r.readline();
	data += line;

	t = header.split("#");
	tt = t[1].split("/");
	barcode = tt[0]
	o.write(t[0] + "/" + tt[1].split("\t")[0].split(" ")[0] + "\tBX:Z:" + barcode + "\n");
	o.write(data);

	header = r.readline()[:-1]

r.close();
o.close();
