#!/usr/bin/python

import os
	
# multimap (dictionary with lists of values for each key)
class KVMap:

	def __init__(self,str=None):
		
		if isinstance(str,KVMap):
			self.dat = dict(str.dat)
		else:
			self.dat = {}
			if str is not None:
				for p in [w.split('=') for w in str.split('\t') ]:
					if len(p) != 2:
						continue
					self.insert(p[0].strip(),p[1].strip())
	
	# insert new entry for given key
	def insert(self,key,value):
		self.dat[key] = self.dat.get(key,[]) + [value,]
	
	# get first value for given key
	def getFirst(self,k,default=None):
		v = self.dat.get(k,[])
		if not len(v):
			return default
		return v[0]
	
	# get first value for given key as a float
	def getFirstF(self,k,default=None):
		if default is not None:
			return float(self.getFirst(k,str(default)))
		else:
			s = self.getFirst(k,None)
			if s is None:
				return None
			return float(s)
			
	# set attributes for float values
	def loadFloats(self,names):
		for nm in names:
			self.__dict__[nm] = self.getFirstF(nm)
	
	# set attributes for string values
	def loadStrings(self,names):
		for nm in names:
			self.__dict__[nm] = self.getFirst(nm)

	# whether this has a matching key:value pair
	def matches(self,k,v):
		return v in self.dat.get(k,[])
	# whether this matches a set of key:value pairs
	def matchesMany(self,kdict):
		for k in kdict:
			if not self.matches(k,kdict[k]):
				return False
		return True
		
	# convert to string
	def toString(self):
		outstr = ""
		for k in self.dat:
			for i in self.dat[k]:
				outstr += str(k)+" = "+str(i)+"\t"
		return outstr[:-1]

# a KVMap string:KVMap
class QFile(KVMap):

	def __init__(self,fname=None):
		self.dat = {}
		self.fname = fname
		
		if fname is not None:
			if not os.path.exists(fname):
				print "No such file",fname
				return
			for l in [z.split(':') for z in open(fname).readlines()]:
				if len(l) < 2:
					continue
				self.dat[l[0]] = self.dat.get(l[0],[]) + [KVMap(l[1]),]
	
	# get first value for key:key
	def getItem(self,k1,k2,default=None):
		m = self.getFirst(k1,None)
		if not m:
			return default
		return m.getFirst(k2,default)
	
	# get first value for key:key as float
	def getItemF(self,k1,k2,default=None):
		m = self.getFirst(k1,None)
		if not m:
			return default
		return m.getFirstF(k2,default)
		
	# get KVMap with matching subkey:value pairs
	def getMatching(self,key,value):
		Q = QFile()
		for k in self.dat.keys():
			for m in self.dat[k]:
				if m.matches(key,value):
					Q.insert(k,m)
		return Q
		
	# get KVMap with matching subkeys:values pairs
	def getMatchingMany(self,requirements):
		Q = QFile()
		for k in self.dat.keys():
			for m in self.dat[k]:
				if m.matchesMany(requirements):
					Q.insert(k,m)
		return Q
	
	# convert to string
	def toString(self):
		outstr = ""
		dkeys = self.dat.keys()
		dkeys.sort()
		for k in dkeys:
			for i in self.dat[k]:
				outstr += str(k)+":\t\t"+i.toString()+"\n"
		return outstr[:-1]
	
