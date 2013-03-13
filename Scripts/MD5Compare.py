#!/usr/bin/python

def md5runfile(f):
	return dict([ [a.split()[1][1:].split("/")[-1], a.split()[0]] for a in open(f).readlines()])
	

def comp_md5(f1,f2):
	d1 = md5runfile(f1)
	d2 = md5runfile(f2)
	d1keys = d1.keys()
	d1keys.sort()

	print "Comparing",len(d1keys),"files\n"
	
	nmatch = 0
	for d in d1keys:
		if d[:3] not in ["run","his"]:
			continue
		#if d[-7:] != ".mid.gz":
		#	continue
		if d[3:5] != "19":
			continue
		if d not in d2:
			print "Missing:",d
		elif d1[d] != d2[d]:
			print "*****",d
		else:
			nmatch += 1

	print "\nMatching: %i"%nmatch

if __name__ == "__main__":
	#comp_md5("/home/mmendenhall/midfiles_2011_LANL_md5.txt","/home/mmendenhall/midfiles_2011_md5.txt")
	comp_md5("/home/mmendenhall/midfiles_2011_LANL_md5.txt","/home/ucna/md5check_sum_ncsu_raw.txt")