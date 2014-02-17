#!/usr/bin/python

#import library needed for bash stuff
import os
import os.path

newoctet = ""
i = 0
#convert runlog into octet list format
with open("OctetList_20112012.txt", "wt") as out:
	for line in open("UCNA Run Log.txt"):
		if (line[:4] == '### '):
			if (i>0):
				out.write(newoctet + '\n\n')
			newoctet = 'Octet:\t\t'
			i = 0
		elif (line[:1] == '*' and (line[7:8] == 'B' or line[7:8] == 'A') and (line[9:10] == '0' or line[9:10] == '1' or line[9:10] == '2') and (line[2:3] == '7' or line[2:3] == '8' or line[2:3] == '9') and line[8:9] != 'g'):
                        newoctet += line[7:10] + ' = ' + line[1:6] + '\t'
                        i += 1
		elif (line[:1] == '*' and (line[7:8] == 'B' or line[7:8] == 'A') and (line[2:3] == '7' or line[2:3] == '8' or line[2:3] == '9') and line[8:9] != 'g'):
			newoctet += line[7:9] + ' = ' + line[1:6] + '\t'
			i += 1
	out.write(newoctet + '\n\n')
