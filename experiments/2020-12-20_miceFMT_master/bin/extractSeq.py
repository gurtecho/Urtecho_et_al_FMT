#!/usr/bin/env python
import os
import sys
def main():
	In = sys.argv[1]
	ref = sys.argv[2]
	Out = sys.argv[3]
	f = open(ref,"rb")
	data = f.readlines()
	f.close()
	Pool = {}
	tmpSeq = ""
	tmpLabel = ""
	flag = 0
	for each in data:
		each = each[:-1]
		if each[0] == ">":
			if flag == 1:
				Pool[tmpLabel] = tmpSeq
			tmpLabel = each[1:]
			tmpSeq = ""
		else:
			tmpSeq += each
		flag = 1
	Pool[tmpLabel] = tmpSeq
	f = open(In,"rb")
	data = f.readlines()
	f.close()
	f = open(Out,"w")
	for each in data:
		each = each[:-1]
		f.writelines([">" + each + os.linesep])
		f.writelines([Pool[each] + os.linesep])
	f.close()
if __name__ == "__main__":
	main()
