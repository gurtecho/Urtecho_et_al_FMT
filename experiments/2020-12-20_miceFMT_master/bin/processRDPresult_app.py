#!/usr/bin/env python
import os
import sys
def main():
	In = sys.argv[1]
	Out = sys.argv[2]
	f = open(In, "rb")
	fw = open(Out, "w")
	tax_list = ["genus", "rootrank", "domain", "phylum", "class", "order", "family"]
	for each in f:
		each = each[:-1]
		tmp = each.split("\t")
		Pool = {}
		for e in tax_list:
			Pool[e] = "NA"
		tmpOtu = tmp[0]
		for i in range((len(tmp) - 1) / 3):
			tmpLabel = tmp[3 * i + 3]
			tmpSeq = tmp[3 * i +  2]
			Pool[tmpLabel] = tmpSeq
		toWrite = [tmpOtu, ] + [Pool[e] for e in tax_list]
		fw.writelines(["\t".join(toWrite) + os.linesep])
	f.close()
	fw.close()
if __name__ == "__main__":
	main()
