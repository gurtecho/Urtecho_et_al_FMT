#!/usr/bin/env python
import os
import sys
def main():
	In = sys.argv[1]
	Out = sys.argv[2]
	f = open(In,"rb")
	data = f.readline()
	f.close()
	Pool = []
	flag = 0
	for e in data:
		if e == "O":
			tmpString = e
			continue
		if e == "t":
			tmpString += e
			continue
		if e == "u":
			tmpString += e
			flag = 1
			continue
		if e != ":":
			if flag == 1:
				tmpString += e
				continue
			else:
				continue
		if e == ":":
			if flag == 1:
				Pool.append(tmpString)
				flag = 0
				continue
			else:
				continue
	f = open(Out,"w")
	for each in Pool:
		f.writelines([each + os.linesep])
if __name__ == "__main__":
	main()
