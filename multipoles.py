import math
import kernels as kl
import terms as tm
import bispec_functions as bf





# l=0

R_0=[tm.R_00]
S_12_0=[tm.S_12_00]
S_23_0=[tm.S_23_00]
S_31_0=[tm.S_31_00]
T_12_0=[tm.T_12_00]
T_23_0=[tm.T_23_00]
T_31_0=[tm.T_31_00]


# l=1

R_1=[]
S_12_1=[]
S_23_1=[]
S_31_1=[]
T_12_1=[]
T_23_1=[]
T_31_1=[]

# l=2

R_2=[tm.R_20,tm.R_21,tm.R_22]
S_12_2=[tm.S_12_20,tm.S_12_21,tm.S_12_22]
S_23_2=[tm.S_23_20,tm.S_23_21,tm.S_23_22]
S_31_2=[tm.S_31_20,tm.S_31_21,tm.S_31_22]
T_12_2=[tm.T_12_20,tm.T_12_21,tm.T_12_22]
T_23_2=[tm.T_23_20,tm.T_23_21,tm.T_23_22]
T_31_2=[tm.T_31_20,tm.T_31_21,tm.T_31_22]

# l=3

R_3=[]
S_12_3=[]
S_23_3=[]
S_31_3=[]
T_12_3=[]
T_23_3=[]
T_31_3=[]

# l=4

R_4=[tm.R_40,tm.R_41,tm.R_42,tm.R_43,tm.R_44]
S_12_4=[tm.S_12_40,tm.S_12_41,tm.S_12_42,tm.S_12_43,tm.S_12_44]
S_23_4=[tm.S_23_40,tm.S_23_41,tm.S_23_42,tm.S_23_43,tm.S_23_44]
S_31_4=[tm.S_31_40,tm.S_31_41,tm.S_31_42,tm.S_31_43,tm.S_31_44]
T_12_4=[tm.T_12_40,tm.T_12_41,tm.T_12_42,tm.T_12_43,tm.T_12_44]
T_23_4=[tm.T_23_40,tm.T_23_41,tm.T_23_42,tm.T_23_43,tm.T_23_44]
T_31_4=[tm.T_31_40,tm.T_31_41,tm.T_31_42,tm.T_31_43,tm.T_31_44]

# l=3

R_5=[]
S_12_5=[]
S_23_5=[]
S_31_5=[]
T_12_5=[]
T_23_5=[]
T_31_5=[]

# l=6

R_6=[tm.R_60,tm.R_61,tm.R_62,tm.R_63,tm.R_64,tm.R_65,tm.R_66]
S_12_6=[tm.S_12_60,tm.S_12_61,tm.S_12_62,tm.S_12_63,tm.S_12_64,tm.S_12_65,tm.S_12_66]
S_23_6=[tm.S_23_60,tm.S_23_61,tm.S_23_62,tm.S_23_63,tm.S_23_64,tm.S_23_65,tm.S_23_66]
S_31_6=[tm.S_31_60,tm.S_31_61,tm.S_31_62,tm.S_31_63,tm.S_31_64,tm.S_31_65,tm.S_31_66]
T_12_6=[tm.T_12_60,tm.T_12_61,tm.T_12_62,tm.T_12_63,tm.T_12_64,tm.T_12_65,tm.T_12_66]
T_23_6=[tm.T_23_60,tm.T_23_61,tm.T_23_62,tm.T_23_63,tm.T_23_64,tm.T_23_65,tm.T_23_66]
T_31_6=[tm.T_31_60,tm.T_31_61,tm.T_31_62,tm.T_31_63,tm.T_31_64,tm.T_31_65,tm.T_31_66]

# l=7

R_7=[]
S_12_7=[]
S_23_7=[]
S_31_7=[]
T_12_7=[]
T_23_7=[]
T_31_7=[]

# l=8

R_8=[tm.R_80,tm.R_81,tm.R_82,tm.R_83,tm.R_84,tm.R_85,tm.R_86,tm.R_87,tm.R_88]
S_12_8=[tm.S_12_80,tm.S_12_81,tm.S_12_82,tm.S_12_83,tm.S_12_84,tm.S_12_85,tm.S_12_86,tm.S_12_87,tm.S_12_88]
S_23_8=[tm.S_23_80,tm.S_23_81,tm.S_23_82,tm.S_23_83,tm.S_23_84,tm.S_23_85,tm.S_23_86,tm.S_23_87,tm.S_23_88]
S_31_8=[tm.S_31_80,tm.S_31_81,tm.S_31_82,tm.S_31_83,tm.S_31_84,tm.S_31_85,tm.S_31_86,tm.S_31_87,tm.S_31_88]
T_12_8=[tm.T_12_80,tm.T_12_81,tm.T_12_82,tm.T_12_83,tm.T_12_84,tm.T_12_85,tm.T_12_86,tm.T_12_87,tm.T_12_88]
T_23_8=[tm.T_23_80,tm.T_23_81,tm.T_23_82,tm.T_23_83,tm.T_23_84,tm.T_23_85,tm.T_23_86,tm.T_23_87,tm.T_23_88]
T_31_8=[tm.T_31_80,tm.T_31_81,tm.T_31_82,tm.T_31_83,tm.T_31_84,tm.T_31_85,tm.T_31_86,tm.T_31_87,tm.T_31_88]


# terms ---------------------------------

R=[R_0,R_1,R_2,R_3,R_4,R_5,R_6,R_7,R_8]
S_12=[S_12_0,S_12_1,S_12_2,S_12_3,S_12_4,S_12_5,S_12_6,S_12_7,S_12_8]
S_23=[S_23_0,S_23_1,S_23_2,S_23_3,S_23_4,S_23_5,S_23_6,S_23_7,S_23_8]
S_31=[S_31_0,S_31_1,S_31_2,S_31_3,S_31_4,S_31_5,S_31_6,S_31_7,S_31_8]
T_12=[T_12_0,T_12_1,T_12_2,T_12_3,T_12_4,T_12_5,T_12_6,T_12_7,T_12_8]
T_23=[T_23_0,T_23_1,T_23_2,T_23_3,T_23_4,T_23_5,T_23_6,T_23_7,T_23_8]
T_31=[T_31_0,T_31_1,T_31_2,T_31_3,T_31_4,T_31_5,T_31_6,T_31_7,T_31_8]





# multipoles 



#========================================================
"""Import power spectrum from the bispec_functions.py"""
#========================================================
""""First we import the cosmological parameters"""
#=======================================================
h0=bf.h0
omg_lam0=bf.omg_lam0
omg_m0=bf.omg_m0
omg_k0=bf.omg_k0
#=======================================================
"""Tracer Power spectrum"""
#=======================================================
"""First set a value of redshift z."""
zz=0.0
"""The growing mode, growth factor and CDM power-spectrum is"""
Dgrow, flog, Pk=bf.initialize(h0,omg_lam0,omg_m0,omg_k0,zz)
#=======================================================
def Pt(kk,b):
	"""This returns the tracer power-spectrum is"""
	return b*b*Pk(kk)
#========================================================
"""Real space bi-spectrum"""
#========================================================


def bispec_real(k1,mu,tt,b,gamma):
	k2=tt*k1
	k3=k1*kl.sfac(mu,tt)
	pk1=Pt(k1,b)
	pk2=Pt(k2,b)
	pk3=Pt(k3,b)
	bs=(2./b)*((kl.F12(mu,tt) + gamma/2.)*pk1*pk2 + (kl.F23(mu,tt) + gamma/2.)*pk2*pk3 + (kl.F31(mu,tt) + gamma/2.)*pk3*pk1)
	return bs

#========================================================
"""Real space bi-spectrum"""
#========================================================


def multipole(l,m,k1,mu,tt,betav,b,gamma):
        k2=tt*k1
        k3=k1*kl.sfac(mu,tt)
        pk1=Pt(k1,b)
        pk2=Pt(k2,b)
        pk3=Pt(k3,b)
        beta=betav
        if betav<0:
                beta=flog/b
                print('beta=',beta)

        G12=kl.G12(mu,tt)
        G23=kl.G23(mu,tt)
        G31=kl.G31(mu,tt)
        dG12=kl.dG12(mu,tt)
        dG23=kl.dG23(mu,tt)
        dG31=kl.dG31(mu,tt)
        Rv = R[l][m](mu,tt,beta)
        S12v= S_12[l][m](mu,tt,beta)
        S23v= S_23[l][m](mu,tt,beta)
        S31v= S_31[l][m](mu,tt,beta)
        T12v= T_12[l][m](mu,tt,beta)
        T23v= T_23[l][m](mu,tt,beta)
        T31v= T_31[l][m](mu,tt,beta)
        bs=(2./b)*(Rv*G12 + S12v*(dG12+gamma/2.)-b*T12v)*pk1*pk2
        bs+=(2./b)*(Rv*G23 + S23v*(dG23+gamma/2.)-b*T23v)*pk2*pk3
        bs+=(2./b)*(Rv*G31 + S31v*(dG31+gamma/2.)-b*T31v)*pk3*pk1
        return bs
