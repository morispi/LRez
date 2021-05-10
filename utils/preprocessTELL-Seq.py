import sys
import re
import csv
import os

r1 = open(sys.argv[1]);
r2 = open(sys.argv[2]);
i1 = open(sys.argv[3]);

o1 = open(os.path.splitext(sys.argv[1])[0] + "_barcodes_in_headers" + os.path.splitext(sys.argv[1])[1], "w");
o2 = open(os.path.splitext(sys.argv[2])[0] + "_barcodes_in_headers" + os.path.splitext(sys.argv[2])[1], "w");

header1 = r1.readline()[:-1];
while header1 != '':
	line = r1.readline();
	data1 = line;
	line = r1.readline();
	data1 += line;
	line = r1.readline();
	data1 += line;

	header2 = r2.readline()[:-1];
	line = r2.readline();
	data2 = line;
	line = r2.readline();
	data2 += line;
	line = r2.readline();
	data2 += line;

	header3 = i1.readline();
	barcode = i1.readline();
	i1.readline();
	i1.readline();

	h = header1.split(" ")[0].split("\t")[0];
	o1.write(h + "\tBX:Z:" + barcode);
	o1.write(data1);
	h = header2.split(" ")[0].split("\t")[0];
	o2.write(h + "\tBX:Z:" + barcode);
	o2.write(data2);

	header1 = r1.readline()[:-1]

r1.close();
r2.close();
i1.close();
o1.close();
o2.close();
