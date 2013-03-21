from ROOT import TTree, TFile

import sys

fname = sys.argv[1]
f = open(fname)

header = []
for line in f:
	if line.strip().find('#') >= 0: 
		header = line.strip()[1:].split('\t')
	else:
		break

tfile = TFile(sys.argv[2], "create")
ttree = TTree(sys.argv[3], "")

from array import array
s = array("f", [0]*len(header))
header[0] += "/F"
ttree.Branch("vals", s, ":".join(header))

while True:
	for i,v in enumerate(line.strip().split('\t')):
		s[i] = float(v)
		
	ttree.Fill()
	print ".",
	try:
		line = f.next()
	except:
		break

ttree.Write()
tfile.Close()
