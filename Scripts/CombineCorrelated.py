#!/usr/bin/python

# implementation of eLog 396 correlated error combination formalism

from math import *
from numpy import zeros,matrix

class MeasCombiner:
	
	def __init__(self):
		self.mu = []
		self.err = []
		self.corrs = {}
	
	def new_meas_mu(self,m):
		self.mu.append(m)
		self.err.append([])
		return len(self.mu)-1

	def new_meas_err(self,meas_id,err):
		self.err[meas_id].append(err)
		err_id = (meas_id,len(self.err[meas_id])-1)
		self.corrs[(err_id,err_id)] = 1.0
		return err_id

	def add_correlation(self,err_id_1,err_id_2,corr):
		self.corrs[(err_id_1,err_id_2)] = corr
		self.corrs[(err_id_2,err_id_1)] = corr

	def calc_combo(self):
		print "Corrs = ",self.corrs
		
		self.ep_i = []
		n = len(self.mu)
		for i in range(n):
			self.ep_i.append(0)
			for j in range(len(self.err[i])):
				for k in range(len(self.err[i])):
					self.ep_i[-1] += self.err[i][j]*self.err[i][k]*self.corrs.get(((i,j),(i,k)),0)
			self.ep_i[-1] = sqrt(self.ep_i[-1])
		print "ep_i =",self.ep_i

		self.C_ij = matrix(zeros((n,n)))
		self.M_ij = matrix(zeros((n,n)))
		for i in range(n):
			for j in range(n):
				for k in range(len(self.err[i])):
					for l in range(len(self.err[j])):
						self.M_ij[i,j] += self.err[i][k]*self.err[j][l]*self.corrs.get(((i,k),(j,l)),0)
				self.C_ij[i,j] = self.M_ij[i,j] / (self.ep_i[i]*self.ep_i[j])
		print "C_ij =",self.C_ij
		print "M_ij =",self.M_ij

		self.W = self.M_ij.getI() * matrix([[1] for i in range(n)])
		self.W /= sum(self.W)
		print "weight =",self.W

		err = 0
		combo = 0
		for i in range(n):
			combo += self.mu[i]*self.W[i,0]
			for j in range(n):
				err += self.W[i,0]*self.ep_i[i]*self.C_ij[i,j]*self.ep_i[j]*self.W[j,0]
		err = sqrt(err)
		print
		print "err =",err
		print "combo =",combo

	def errcombo(self,errlist):
		err = 0
		n = len(self.mu)
		for i in range(n):
			for j in range(n):
				for k in range(len(self.err[i])):
					for l in range(len(self.err[j])):
						if (i,k) in errlist and (j,l) in errlist:
							err += self.W[i,0]*self.err[i][k]*self.W[j,0]*self.err[j][l]*self.corrs.get(((i,k),(j,l)),0)
		err = sqrt(err)
		return err

if __name__=="__main__":
	# eLog 396 example
	if 0:
		MC = MeasCombiner()
		mu0 = MC.new_meas_mu(100.0)
		MC.new_meas_err(mu0,0.6)
		err01 = MC.new_meas_err(mu0,0.3)
		err02 = MC.new_meas_err(mu0,0.2)
		mu1 = MC.new_meas_mu(100.0)
		MC.new_meas_err(mu1,0.6)
		err11 = MC.new_meas_err(mu1,0.3)
		err12 = MC.new_meas_err(mu1,0.6)
		MC.add_correlation(err01,err11,1.0)
		MC.add_correlation(err02,err12,1.0)
		MC.calc_combo()

	# 2012 Before + After
	if 1:
		MC = MeasCombiner()
		mu0 = MC.new_meas_mu(-0.11851)
		stat0 = MC.new_meas_err(mu0,0.00075)	# stat
		err01 = MC.new_meas_err(mu0,0.000712) # syst
		dst0 = MC.new_meas_err(mu0,0.000119) # depol stat
		err02 = MC.new_meas_err(mu0,0.000640) # depol syst
		
		mu1 = MC.new_meas_mu(-0.12100)
		stat1 = MC.new_meas_err(mu1,0.00078) # stat
		err11 = MC.new_meas_err(mu1,0.000726) # syst
		dst1 = MC.new_meas_err(mu1,0.000387) # depol stat
		err12 = MC.new_meas_err(mu1,0.000666) # depol syst
		
		MC.add_correlation(err01,err11,1.0)
		MC.add_correlation(err02,err12,1.0)
		
		MC.calc_combo()
		print "stat =",MC.errcombo([stat0,stat1])
		print "syst =",MC.errcombo([err01,err11])
		print "depol =",MC.errcombo([dst0,dst1,err02,err12])
		print "systot =",MC.errcombo([err01,err11,dst0,dst1,err02,err12])
		print
	
	# 2010 + 2012 results, + side
	if 1:
		MC = MeasCombiner()
		mu0 = MC.new_meas_mu(0.11954)
		MC.new_meas_err(mu0,0.00055)	# stat
		err01 = MC.new_meas_err(mu0,0.00072) # syst
		err02 = MC.new_meas_err(mu0,0.00067) # depol
		mu1 = MC.new_meas_mu(0.11966)
		MC.new_meas_err(mu1,0.00089) # stat
		err11 = MC.new_meas_err(mu1,0.001254) # syst
		err12 = MC.new_meas_err(mu1,0.000622) # depol
		MC.add_correlation(err01,err11,1.0)
		MC.calc_combo()
		print

	# 2010 + 2012 results, - side
	if 1:
		MC = MeasCombiner()
		mu0 = MC.new_meas_mu(0.11954)
		MC.new_meas_err(mu0,0.00055)	# stat
		err01 = MC.new_meas_err(mu0,0.00072) # syst
		err02 = MC.new_meas_err(mu0,0.00067) # depol
		mu1 = MC.new_meas_mu(0.11966)
		MC.new_meas_err(mu1,0.00089) # stat
		err11 = MC.new_meas_err(mu1,0.00123) # syst
		MC.add_correlation(err01,err11,1.0)
		MC.calc_combo()
		print
