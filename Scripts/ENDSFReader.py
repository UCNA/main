#!/usr/bin/python

def parse_n_uncert(d,nuncdig=2):
	nstr = d[:-nuncdig].strip().upper()
	if not nstr.strip():
		return None
	n = float(nstr)
	expo = 0
	if 'E' in nstr:
		expo = int(nstr.split('E')[1])
		nstr = nstr.split('E')[0]
	if '.' in nstr:
		expo -= len(nstr.split('.'))
	uncert = 0
	if d[-nuncdig:].strip():
		uncert = int(d[-nuncdig:].strip())*(10**expo)
	return (n,uncert)
	

class ENDSFLine:
	def __init__(self,l):
		self.dat = l
		self.NUCID = l[:5].strip()
		self.F6 = l[5]
		self.F7 = l[6]
		self.F8 = l[7]
		self.F9 = l[8]
		self.conts = l[9:-1]

	def __repr__(self):
		return '<[%s]%s%s%s: %s>'%(self.NUCID,self.F6,self.F7,self.F8,self.conts)

class ENDSF_Norm_Card:
	def __init__(self,l):
		self.l = l
		self.NR = parse_n_uncert(l.dat[9:21])		# conversion from RI in gamma record to photons per 100 neutron captures
		self.NT = parse_n_uncert(l.dat[21:31])		# conversion from TI in gamma record to transitions per 100 parent decays/captures
		self.BR = parse_n_uncert(l.dat[31:42])		# branching ratio factor from intensity per 100 decays in this branch to intensity per 100 decays of parent
		self.NB = parse_n_uncert(l.dat[41:55],6)	# conversion from IB in B-, IB,IE,TI in EC record to intensities per 100 decays in this branch
		self.NP = parse_n_uncert(l.dat[55:64])		# conversion from 100 delayed transition intensities to 100 precursor decays

class ENDSF_Norm_Prod_Card:
	def __init__(self,l):
		self.l = l
		self.NRxBR = parse_n_uncert(l.dat[9:21])	# conversion from RI in gamma record to photons per 100 parent decays
		self.NTxBR = parse_n_uncert(l.dat[21:31])	# conversion from TI in gamma record to transitions per 100 parent decays
		self.NBxBR = parse_n_uncert(l.dat[41:55],6)	# conversion from IB in B-, IB,IE,TI in EC record to intensities per 100 parent decays
		self.NP = parse_n_uncert(l.dat[55:62],0)	# conversion from 100 delayed transition intensities to 100 precursor decays
		self.COM = l.dat[76].strip()
		self.OPT = l.dat[77].strip()
		

	def __repr__(self):
		s = "<ProdNorm:"
		if self.NRxBR:
			s += " NRxBR=%g~%g"%self.NRxBR
		if self.NTxBR:
			s += " NTxBR=%g~%g"%self.NTxBR
		if self.NBxBR:
			s += " NBxBR=%g~%g"%self.NBxBR
		if self.NP:
			s += " NP=%g~%g"%self.NP
		if self.OPT:
			s += " OPT="+self.OPT
		return s+">"


class ENDSF_Level_Card:
	def __init__(self,l):
		self.l = l
		self.E = parse_n_uncert(l.dat[9:21])
		self.J = l.dat[21:39].strip()
		self.T = l.dat[39:55].strip() # TODO parse time units
		self.L = l.dat[55:64].strip()
		self.S = parse_n_uncert(l.dat[64:76])
		self.C = l.dat[76]
		self.MS = l.dat[77:79].strip()
		self.Q = l.dat[79].strip()
			
	def __repr__(self):
		return "<Lvl: %g~%g %s%s>"%(self.E[0],self.E[1],self.J,self.Q)

class ENDSF_Gamma_Card:
	def __init__(self,l):
		self.l = l
		self.E = parse_n_uncert(l.dat[9:21])		# gamma energy in keV
		self.RI = parse_n_uncert(l.dat[21:31])		# relative photon intensity
		self.M = l.dat[31:41].strip()				# multipolarity
		self.MR = parse_n_uncert(l.dat[41:55],6)	# mixing ratio
		self.CC = parse_n_uncert(l.dat[55:64])		# total conversion coefficient
		self.TI = parse_n_uncert(l.dat[64:76])		# relative total transition intensity
		self.C = l.dat[76].strip()					# comment flag
		self.COIN = l.dat[77].strip()				# confirmed by coincidence?
		self.Q = l.dat[79].strip()					# placement uncertainty
		self.level_from = None
	
	def __repr__(self):
		s = "<G(%i)"%self.level_from + " %g~%g"%self.E
		if self.RI:
			s += " RI=%g~%g"%self.RI
		if self.TI:
			s += " TI=%g~%g"%self.TI
		if self.M:
			s += " M=%s"%self.M
		if self.C:
			s += " C=%s"%self.C
		if self.Q:
			s += " Q=%s"%self.Q
		return s+">"

class ENDSFReader:
	def __init__(self,fname):
		self.data = [ENDSFLine(l) for l in open(fname,"r").readlines() if len(l.strip())]
		
		self.ID_rec = self.data[0]
		self.Hist_recs = []
		self.Comment_recs = []
		self.Parent_recs = []
		self.Norm_recs = []
		self.Level_recs = []
		self.Gamma_recs = []
		
		n = 1
		while n<len(self.data):
			
			# continuation records
			if self.data[n].F6 != ' ':
				print "Unhandled continuation:",self.data[n]
			
			elif self.data[n].F8 == 'H':
				if self.data[n].F6 == ' ':
					self.Hist_recs.append(self.data[n])
				else:
					self.Hist_recs[-1].conts += self.data[n].conts
		
			elif self.data[n].F7 in ['c','C','d','D','t','T']:
				self.Comment_recs.append(self.data[n])
			
			elif self.data[n].F8 == 'P':
				self.Parent_recs.append(self.data[n])
										
			elif self.data[n].F8 == 'N' and self.data[n].F7 in " P":
				if self.data[n].F7 == ' ':
					self.Norm_recs.append(ENDSF_Norm_Card(self.data[n]))
				elif self.data[n].F7 == 'P':
					self.Norm_recs.append(ENDSF_Norm_Prod_Card(self.data[n]))
						
			elif self.data[n].F8 == 'L':
				self.Level_recs.append(ENDSF_Level_Card(self.data[n]))

			elif self.data[n].F8 == 'G':
				self.Gamma_recs.append(ENDSF_Gamma_Card(self.data[n]))
				self.Gamma_recs[-1].level_from = len(self.Level_recs)-1
						
			else:
				print "Unknown card:",self.data[n]
						
			n += 1
		
			


if __name__ == "__main__":
	E = ENDSFReader("/Users/michael/Desktop/ENDSF_Cu65.txt")
	print E.Norm_recs
	for g in E.Gamma_recs:
		if g.RI and g.RI[0]>0.1:
			print g
