from __future__ import unicode_literals
import numpy as np
import math
from scipy.interpolate import griddata
from scipy.interpolate import interp1d



#****************** COSMOLOGY CALCULATIONS *****************#
#----------- read cosmological parameters ----------------
h0=0.67
H0=h0*100.
omg_m0=0.31
omg_lam0=0.69
omg_k0=(1.-omg_m0-omg_lam0)
#---------------------------------------------------------

#----------- Hubble parameter ----------------
def Ea(af,omgL,omgM,omgK):
	return (omgM*pow(af,-3.)+omgK*pow(af,-2.)+omgL)**0.5


def Ha(af,omgL,omgM,omgK):
	return H0*Ea(af,omgL,omgM,omgK)

#--------------- Deceleration parameter -----------------
def Qa(af,omgL,omgM,omgK):
	return (omgM*pow(af,-3.)-2.*omgL)/(2.*Ea(af,omgL,omgM,omgK)**2.)

#------------- Growth factor ---------------------------
def IDf(af,omgL,omgM,omgK):
	return pow(af*Ea(af,omgL,omgM,omgK),-3.)

def Da(af,omgL,omgM,omgK):
	a=0.00001
	b=af
	n=1000

	h = float(b - a)/n

	s = 0.0
	s += IDf(a,omgL,omgM,omgK)/2.0

	for i in range(1, n):
		s += IDf(a+i*h,omgL,omgM,omgK)

	s += IDf(b,omgL,omgM,omgK)/2.0
	
	return Ea(af,omgL,omgM,omgK)*s*h

#------------------- Logarithmic growth rate ------------
def fa(af,omgL,omgM,omgK):
	return 1./(af*af*Ea(af,omgL,omgM,omgK)*Ea(af,omgL,omgM,omgK)*Da(af,omgL,omgM,omgK)) - (1.+Qa(af,omgL,omgM,omgK))

#------------------- Initialize -------------------------
def initialize(hh,omgL,omgM,omgK,zz):
	scalef=1./(1.+zz)
	Dz=Da(scalef,omgL,omgM,omgK)
	fz=fa(scalef,omgL,omgM,omgK)

	Pr=np.loadtxt("./pk.dat")
	for i in range (len(Pr)):
		Pr[i,0]=Pr[i,0]*hh
		Pr[i,1]=Pr[i,1]*Dz*Dz/(hh**3.)

	Pz = interp1d(Pr[:,0], Pr[:,1], kind = 'cubic')
	return Dz, fz, Pz





