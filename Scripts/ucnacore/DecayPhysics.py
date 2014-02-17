#!/usr/bin/python

# various constants and utilities for neutron decay

from math import *

neutronBetaEp = 782.347				# neutron beta decay endpoint, keV
m_e = 511.00						# electron mass, keV/c^2
m_p = 938272.046					# proton mass, keV/c^2
lmbda = abs(-1.2694) 				# +/-0.0028, PDG 2010 value lambda, Wilkinson sign convention
A0_PDG = -0.1173					# +/-0.0013, PDG 2010 value

def beta(KE, m = m_e):
	"""beta = v/c as a function of kinetic energy"""
	return sqrt(KE**2+2*m_e*KE)/(m_e+KE)

def p_e(KE):
	"""electron momentum as a function of kinetic energy"""
	return sqrt(KE**2 + 2*m_e*KE)

def S0(KE):
	"""beta spectrum phase space as a function of kinetic energy"""
	if not 0 < KE < neutronBetaEp:
		return 0
	return p_e(KE)*(KE+m_e)*(neutronBetaEp-KE)**2

def WilkinsonRWM(KE):
	"""Wilkinson Recoil + Weak Magnetism correction to A"""
	W = (KE+m_e)/m_e
	W0 = neutronBetaEp/m_e+1.
	mu = 2.792847356-(-1.91304273);	# mu_p - mu_n = 2.792847356(23) - -1.91304273(45) PDG 2010
	A_uM = (lmbda+mu)/(lmbda*(1.-lmbda)*(1.+3.*lmbda**2)*m_p/m_e)
	A_1 = lmbda**2 + 2.*lmbda/3. - 1./3.
	A_2 = -lmbda**3 - 3.*lmbda**2 - 5.*lmbda/3. + 1./3.
	A_3 = 2.*lmbda**2*(1.-lmbda)
	return A_uM*(A_1*W0+A_2*W+A_3/W)

def SPol(KE,P):
	"""polarized beta decay spectrum with asymmetry term"""
	return S0(KE)*(1+0.5*P*A0_PDG*beta(KE)*(1+WilkinsonRWM(KE)))


