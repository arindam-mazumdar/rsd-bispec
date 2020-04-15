import math
import kernels as kl
Power=pow
Sqrt=math.sqrt
Pi=math.pi



#==Starting l=0,m=0


def R_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (7*(15 + beta*(10 + beta + 2*beta*Power(mu,2)))*Power(ss,2) + beta*(-2*(35 + 42*beta + 9*Power(beta,2))*mu*tt - 12*Power(beta,2)*Power(mu,3)*tt + (35 + beta*(28 + 3*beta))*(1 + Power(tt,2)) + 2*beta*(7 + 6*beta)*Power(mu,2)*(1 + Power(tt,2))))/(105.*Power(ss,2))

def S_12_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (15 + beta*(10 + beta + 2*beta*Power(mu,2)))/15.

def S_23_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (5*(3 + beta)*Power(ss,2) + beta*(5 + beta + 2*beta*Power(mu,2) - 2*(5 + 3*beta)*mu*tt + (5 + 3*beta)*Power(tt,2)))/(15.*Power(ss,2))

def S_31_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (5*(3 + beta)*Power(ss,2) + beta*(5 + 3*beta - 2*(5 + 3*beta)*mu*tt + (5 + beta + 2*beta*Power(mu,2))*Power(tt,2)))/(15.*Power(ss,2))

def T_12_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(-3*(35 + beta*(35 + beta*(9 + beta))) - 12*beta*(1 + beta)*(7 + 2*beta)*Power(mu,2) - 8*Power(beta,3)*Power(mu,4) + 3*(1 + beta)*(35 + beta*(28 + 5*beta))*mu*tt + 4*Power(beta,2)*(9 + 5*beta)*Power(mu,3)*tt))/630. + (beta*(3*(1 + beta)*(35 + beta*(28 + 5*beta))*mu + 4*Power(beta,2)*(9 + 5*beta)*Power(mu,3) + (-3*(35 + beta*(35 + beta*(9 + beta))) - 12*beta*(1 + beta)*(7 + 2*beta)*Power(mu,2) - 8*Power(beta,3)*Power(mu,4))*tt))/(630.*tt)

def T_23_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(21*(5 + 3*beta)*mu*Power(ss,4) + 6*beta*Power(ss,2)*(3*mu*(7 + 3*beta + 2*beta*Power(mu,2)) - 2*(7 + 3*beta + 14*Power(mu,2) + 12*beta*Power(mu,2))*tt + 3*(7 + 5*beta)*mu*Power(tt,2)) + Power(beta,2)*(5*mu*(9 + beta*(3 + 4*Power(mu,2))) - 4*(9 + 36*Power(mu,2) + beta*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*tt + 6*mu*(9*(3 + 2*Power(mu,2)) + 5*beta*(3 + 4*Power(mu,2)))*Power(tt,2) - 4*(9 + 36*Power(mu,2) + 5*beta*(1 + 6*Power(mu,2)))*Power(tt,3) + 5*(9 + 7*beta)*mu*Power(tt,4))))/(630.*Power(ss,4)*tt) + (beta*(beta*(-63 - 18*beta*(1 + 4*Power(mu,2)) - Power(beta,2)*(3 + 24*Power(mu,2) + 8*Power(mu,4)) + 3*mu*(63 + 18*beta*(3 + 2*Power(mu,2)) + 5*Power(beta,2)*(3 + 4*Power(mu,2)))*tt - 3*(21 + 42*Power(mu,2) + 18*beta*(1 + 4*Power(mu,2)) + 5*Power(beta,2)*(1 + 6*Power(mu,2)))*Power(tt,2) + (63 + 90*beta + 35*Power(beta,2))*mu*Power(tt,3)) + 3*Power(ss,2)*(35*(-1 + mu*tt) + 14*beta*(-1 - 2*Power(mu,2) + 3*mu*tt) + 3*Power(beta,2)*(-1 - 4*Power(mu,2) + 5*mu*tt))))/(630.*Power(ss,4))

def T_31_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*tt*(3*Power(ss,2)*(-((35 + 42*beta + 15*Power(beta,2))*mu) + (35 + 14*beta*(1 + 2*Power(mu,2)) + 3*Power(beta,2)*(1 + 4*Power(mu,2)))*tt) + beta*(-((63 + 90*beta + 35*Power(beta,2))*mu) + 3*(21 + 42*Power(mu,2) + 18*beta*(1 + 4*Power(mu,2)) + 5*Power(beta,2)*(1 + 6*Power(mu,2)))*tt - 3*mu*(63 + 18*beta*(3 + 2*Power(mu,2)) + 5*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,2) + (63 + 18*beta*(1 + 4*Power(mu,2)) + Power(beta,2)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,3))))/(630.*Power(ss,4)) - (beta*tt*(21*(5 + 3*beta)*mu*Power(ss,4) + 6*beta*Power(ss,2)*(3*(7 + 5*beta)*mu - 2*(7 + 3*beta + 14*Power(mu,2) + 12*beta*Power(mu,2))*tt + 3*mu*(7 + 3*beta + 2*beta*Power(mu,2))*Power(tt,2)) + Power(beta,2)*(5*(9 + 7*beta)*mu - 4*(9 + 36*Power(mu,2) + 5*beta*(1 + 6*Power(mu,2)))*tt + 6*mu*(9*(3 + 2*Power(mu,2)) + 5*beta*(3 + 4*Power(mu,2)))*Power(tt,2) - 4*(9 + 36*Power(mu,2) + beta*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,3) + 5*mu*(9 + beta*(3 + 4*Power(mu,2)))*Power(tt,4))))/(630.*Power(ss,4))

#==End l=0,m=0




#==Starting l=2,m=0


def R_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(14 + beta*(13 + beta + (11 + 9*beta)*Power(mu,2)) + (7 + beta + (21 + 11*beta)*Power(mu,2))*Power(ss,2) - 2*mu*(14 + 3*beta*(5 + beta) + beta*(9 + 7*beta)*Power(mu,2))*tt + (-7 - 5*beta + (21 + beta*(29 + 6*beta))*Power(mu,2) + 4*Power(beta,2)*Power(mu,4))*Power(tt,2)))/(21.*Power(ss,2))

def S_12_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(7 + beta + 21*Power(mu,2) + 11*beta*Power(mu,2)))/21.

def S_23_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(14 + beta + 11*beta*Power(mu,2) + 7*(-1 + 3*Power(mu,2))*Power(ss,2) - 2*mu*(14 + 3*beta + 9*beta*Power(mu,2))*tt + (7 + 6*beta)*(-1 + 3*Power(mu,2))*Power(tt,2)))/(21.*Power(ss,2))

def S_31_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(2*(7 + 6*beta) + 14*Power(ss,2) - 4*(7 + 6*beta)*mu*tt + (-7 + beta + 21*Power(mu,2) + 11*beta*Power(mu,2))*Power(tt,2)))/(21.*Power(ss,2))

def T_12_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-3*(154 + beta*(297 + beta*(176 + 25*beta)))*mu - beta*(297 + beta*(462 + 205*beta))*Power(mu,3) + (3*(-77 + beta*(-44 + beta*(11 + 2*beta))) + 3*(231 + beta*(440 + beta*(231 + 46*beta)))*Power(mu,2) + 8*Power(beta,2)*(33 + 17*beta)*Power(mu,4))*tt))/(1386.*tt) + (beta*(231*(-1 + mu*tt) + 33*beta*(-7 - 11*Power(mu,2) + 9*(mu + Power(mu,3))*tt) + 33*Power(beta,2)*(-1 - 12*Power(mu,2) - 2*Power(mu,4) + 3*(mu + 4*Power(mu,3))*tt) + Power(beta,3)*(-3 - 69*Power(mu,2) - 68*Power(mu,4) + 5*mu*(3 + 19*Power(mu,2) + 6*Power(mu,4))*tt)))/693.

def T_23_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(33*mu*(14 + beta*(3 + 9*Power(mu,2)))*Power(ss,4) + 66*beta*Power(ss,2)*(mu*(12 + beta*(3 + 7*Power(mu,2))) - 2*(1 + (11 + 6*beta)*Power(mu,2) + 4*beta*Power(mu,4))*tt + mu*(3 + (9 + 10*beta)*Power(mu,2))*Power(tt,2)) + Power(beta,2)*(5*mu*(66 + beta*(15 + 41*Power(mu,2))) - 4*(33*(1 + 9*Power(mu,2)) + 2*beta*(3 + 69*Power(mu,2) + 68*Power(mu,4)))*tt + 6*mu*(33*(3 + 7*Power(mu,2)) + 10*beta*(3 + 19*Power(mu,2) + 6*Power(mu,4)))*Power(tt,2) - 4*(66*Power(mu,2)*(3 + 2*Power(mu,2)) + 5*beta*(-1 + 21*Power(mu,2) + 36*Power(mu,4)))*Power(tt,3) + 5*(66*Power(mu,3) + 7*beta*mu*(-1 + 9*Power(mu,2)))*Power(tt,4))))/(1386.*Power(ss,4)*tt) + (beta*(beta*(-2*(198 + 33*beta*(1 + 9*Power(mu,2)) + Power(beta,2)*(3 + 69*Power(mu,2) + 68*Power(mu,4))) + 6*mu*(198 + 33*beta*(3 + 7*Power(mu,2)) + 5*Power(beta,2)*(3 + 19*Power(mu,2) + 6*Power(mu,4)))*tt - 3*(33 + 363*Power(mu,2) + 132*beta*Power(mu,2)*(3 + 2*Power(mu,2)) + 5*Power(beta,2)*(-1 + 21*Power(mu,2) + 36*Power(mu,4)))*Power(tt,2) + mu*(660*beta*Power(mu,2) + 99*(1 + 3*Power(mu,2)) + 35*Power(beta,2)*(-1 + 9*Power(mu,2)))*Power(tt,3)) + 66*Power(ss,2)*(-7 + 7*mu*tt + Power(beta,2)*Power(mu,2)*(-3 - 2*Power(mu,2) + 5*mu*tt) + beta*(-1 - 11*Power(mu,2) + 3*mu*tt + 9*Power(mu,3)*tt))))/(1386.*Power(ss,4))

def T_31_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*tt*(33*Power(ss,2)*(-2*(7 + 12*beta + 5*Power(beta,2))*mu + (-7 + 2*beta + Power(beta,2) + 21*Power(mu,2) + 22*beta*Power(mu,2) + 9*Power(beta,2)*Power(mu,2))*tt) + beta*(-4*(99 + 165*beta + 70*Power(beta,2))*mu + 3*(33 + 363*Power(mu,2) + 66*beta*(1 + 9*Power(mu,2)) + 5*Power(beta,2)*(5 + 51*Power(mu,2)))*tt - 3*mu*(99*(1 + 3*Power(mu,2)) + 66*beta*(3 + 7*Power(mu,2)) + 5*Power(beta,2)*(15 + 41*Power(mu,2)))*Power(tt,2) + 2*(66*beta*Power(mu,2)*(3 + 2*Power(mu,2)) + 99*(-1 + 3*Power(mu,2)) + Power(beta,2)*(3 + 69*Power(mu,2) + 68*Power(mu,4)))*Power(tt,3))))/(1386.*Power(ss,4)) - (beta*tt*(33*(7 + 6*beta)*mu*Power(ss,4) + 33*beta*Power(ss,2)*(2*(6 + 5*beta)*mu - 2*(1 + beta + 11*Power(mu,2) + 9*beta*Power(mu,2))*tt + mu*(3 + 3*beta + 9*Power(mu,2) + 7*beta*Power(mu,2))*Power(tt,2)) + Power(beta,2)*(5*(33 + 28*beta)*mu - 2*(33*(1 + 9*Power(mu,2)) + 5*beta*(5 + 51*Power(mu,2)))*tt + 3*mu*(99 + 75*beta + 231*Power(mu,2) + 205*beta*Power(mu,2))*Power(tt,2) - 4*(33*Power(mu,2)*(3 + 2*Power(mu,2)) + beta*(3 + 69*Power(mu,2) + 68*Power(mu,4)))*Power(tt,3) + 5*(33*Power(mu,3) + beta*mu*(3 + 19*Power(mu,2) + 6*Power(mu,4)))*Power(tt,4))))/(693.*Power(ss,4))

#==End l=2,m=0




#==Starting l=2,m=1


def R_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.6666666666666666)*beta*Sqrt(1 - Power(mu,2))*(mu*(beta*(9 + 5*beta) + 3*(7 + 3*beta)*Power(ss,2)) - 3*(7 + beta*(6 + beta + 2*(3 + 2*beta)*Power(mu,2)))*tt + mu*(21 + beta*(27 + beta*(6 + 4*Power(mu,2))))*Power(tt,2)))/(21.*Power(ss,2))

def S_12_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.6666666666666666)*beta*(7 + 3*beta)*mu*Sqrt(1 - Power(mu,2)))/7.

def S_23_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.6666666666666666)*beta*Sqrt(1 - Power(mu,2))*(mu*(3*beta + 7*Power(ss,2)) - (7 + beta*(3 + 6*Power(mu,2)))*tt + (7 + 6*beta)*mu*Power(tt,2)))/(7.*Power(ss,2))

def S_31_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.6666666666666666)*beta*(7 + 3*beta)*Sqrt(1 - Power(mu,2))*tt*(-1 + mu*tt))/(7.*Power(ss,2))

def T_12_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*Sqrt(1 - Power(mu,2))*(-4*beta*mu*(99 + 22*beta*(4 + Power(mu,2)) + 5*Power(beta,2)*(3 + 4*Power(mu,2))) + (231 + 121*Power(beta,2)*(1 + 4*Power(mu,2)) + 99*beta*(3 + 4*Power(mu,2)) + 5*Power(beta,3)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*tt))/(462.*Sqrt(6)) + (beta*Sqrt(1 - Power(mu,2))*(231 - 462*mu*tt + beta*(297 + 121*beta + 15*Power(beta,2) + 6*(33 + beta*(44 + 15*beta))*Power(mu,2) - 2*(396 + beta*(187 + 30*beta))*mu*tt - 16*beta*(11 + 5*beta)*Power(mu,3)*tt)))/(462.*Sqrt(6)*tt)

def T_23_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*Sqrt(1 - Power(mu,2))*(11*Power(ss,2)*(-4*beta*mu*(9 + beta*(3 + 2*Power(mu,2))) + (21 + 18*beta*(1 + 2*Power(mu,2)) + 5*Power(beta,2)*(1 + 4*Power(mu,2)))*tt) + beta*(-20*beta*mu*(11 + beta*(3 + 4*Power(mu,2))) + 3*(99 + 66*beta*(1 + 4*Power(mu,2)) + 5*Power(beta,2)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*tt - 6*mu*(99 + 44*beta*(3 + 2*Power(mu,2)) + 15*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,2) + (99*(1 + 2*Power(mu,2)) + 110*beta*(1 + 4*Power(mu,2)) + 35*Power(beta,2)*(1 + 6*Power(mu,2)))*Power(tt,3))))/(462.*Sqrt(6)*Power(ss,4)) - (beta*Sqrt(1 - Power(mu,2))*(33*(7 + beta*(3 + 6*Power(mu,2)))*Power(ss,4) + 22*beta*Power(ss,2)*(3*(3 + beta + 4*beta*Power(mu,2)) - 4*mu*(9 + 6*beta + 4*beta*Power(mu,2))*tt + (9 + 5*beta + 18*Power(mu,2) + 20*beta*Power(mu,2))*Power(tt,2)) + Power(beta,2)*(55 + 15*beta + 90*beta*Power(mu,2) - 40*mu*(11 + beta*(6 + 8*Power(mu,2)))*tt + 6*(33*(1 + 4*Power(mu,2)) + 5*beta*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,2) - 8*mu*(66 + 44*Power(mu,2) + 15*beta*(3 + 4*Power(mu,2)))*Power(tt,3) + 5*(11 + 44*Power(mu,2) + 7*beta*(1 + 6*Power(mu,2)))*Power(tt,4))))/(462.*Sqrt(6)*Power(ss,4)*tt)

def T_31_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*Sqrt(1 - Power(mu,2))*tt*(11*(21 + 18*beta + 5*Power(beta,2))*Power(ss,2)*(-1 + 2*mu*tt) + beta*(-99 - 110*beta - 35*Power(beta,2) + 6*(99 + 110*beta + 35*Power(beta,2))*mu*tt - 9*(33 + 66*Power(mu,2) + 22*beta*(1 + 4*Power(mu,2)) + 5*Power(beta,2)*(1 + 6*Power(mu,2)))*Power(tt,2) + 4*mu*(99 + 22*beta*(3 + 2*Power(mu,2)) + 5*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,3))))/(462.*Sqrt(6)*Power(ss,4)) - (beta*Sqrt(1 - Power(mu,2))*tt*(33*(7 + 3*beta)*Power(ss,4) + 22*beta*Power(ss,2)*(9 + 5*beta - 4*(9 + 5*beta)*mu*tt + 3*(3 + beta + 6*Power(mu,2) + 4*beta*Power(mu,2))*Power(tt,2)) + Power(beta,2)*(55 + 35*beta - 40*(11 + 7*beta)*mu*tt + 18*(11 + 44*Power(mu,2) + 5*beta*(1 + 6*Power(mu,2)))*Power(tt,2) - 16*mu*(33 + 22*Power(mu,2) + 5*beta*(3 + 4*Power(mu,2)))*Power(tt,3) + 5*(11 + 44*Power(mu,2) + beta*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,4))))/(462.*Sqrt(6)*Power(ss,4))

#==End l=2,m=1




#==Starting l=2,m=2


def R_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2))*(beta*(3 + beta) + 3*(7 + beta)*Power(ss,2) - 6*beta*(3 + beta)*mu*tt + (21 + 21*beta + Power(beta,2)*(2 + 4*Power(mu,2)))*Power(tt,2)))/(21.*Sqrt(6)*Power(ss,2))

def S_12_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(7 + beta)*(-1 + Power(mu,2)))/(7.*Sqrt(6))

def S_23_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2))*(beta + 7*Power(ss,2) - 6*beta*mu*tt + (7 + 6*beta)*Power(tt,2)))/(7.*Sqrt(6)*Power(ss,2))

def S_31_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(7 + beta)*(-1 + Power(mu,2))*Power(tt,2))/(7.*Sqrt(6)*Power(ss,2))

def T_12_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(-1 + Power(mu,2))*(-3*beta*(33 + beta*(22 + 5*beta))*mu + (231 + beta*(264 + beta*(55 + 6*beta + 8*(11 + 3*beta)*Power(mu,2))))*tt))/(462.*Sqrt(6)*tt) - (Power(beta,2)*(-1 + Power(mu,2))*(-33 + 99*mu*tt + beta*(-22 - 3*beta - 2*(11 + 6*beta)*Power(mu,2) + (88 + 15*beta)*mu*tt + 10*beta*Power(mu,3)*tt)))/(231.*Sqrt(6))

def T_23_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-1 + Power(mu,2))*(99*mu*Power(ss,4) + 22*Power(ss,2)*(3*beta*mu - 2*(3 + 2*beta + 4*beta*Power(mu,2))*tt + (9 + 10*beta)*mu*Power(tt,2)) + beta*(15*beta*mu - 4*(11 + 6*beta + 24*beta*Power(mu,2))*tt + 6*mu*(33 + 30*beta + 20*beta*Power(mu,2))*Power(tt,2) - 4*(22 + 15*beta + 44*Power(mu,2) + 60*beta*Power(mu,2))*Power(tt,3) + 5*(22 + 21*beta)*mu*Power(tt,4))))/(462.*Sqrt(6)*Power(ss,4)*tt) - (Power(beta,2)*(-1 + Power(mu,2))*(-2*beta*(11 + 3*beta*(1 + 4*Power(mu,2))) + 6*beta*mu*(33 + 5*beta*(3 + 2*Power(mu,2)))*tt - 3*(33 + 44*beta*(1 + 2*Power(mu,2)) + 15*Power(beta,2)*(1 + 4*Power(mu,2)))*Power(tt,2) + (99 + 220*beta + 105*Power(beta,2))*mu*Power(tt,3) + 22*Power(ss,2)*(-3 + 9*mu*tt + beta*(-1 - 2*Power(mu,2) + 5*mu*tt))))/(462.*Sqrt(6)*Power(ss,4))

def T_31_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(-1 + Power(mu,2))*Power(tt,2)*(11*(21 + beta*(6 + beta))*Power(ss,2) + beta*(99 + 3*beta*(22 + 5*beta) - 9*(33 + beta*(22 + 5*beta))*mu*tt + 2*(99 + beta*(22 + 3*beta + 4*(11 + 3*beta)*Power(mu,2)))*Power(tt,2))))/(462.*Sqrt(6)*Power(ss,4)) + (Power(beta,2)*(-1 + Power(mu,2))*Power(tt,2)*(11*(3 + beta)*Power(ss,2)*(-2 + 3*mu*tt) + beta*(-2*(11 + 5*beta) + 9*(11 + 5*beta)*mu*tt - 4*(11 + 3*beta + 22*Power(mu,2) + 12*beta*Power(mu,2))*Power(tt,2) + 5*mu*(11 + 3*beta + 2*beta*Power(mu,2))*Power(tt,3))))/(231.*Sqrt(6)*Power(ss,4))

#==End l=2,m=2




#==Starting l=4,m=0


def R_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(44 - 16*beta + 4*(33 + 34*beta)*Power(mu,2) + 44*(-1 + 3*Power(mu,2))*Power(ss,2) + 8*mu*(11 + 12*beta - (55 + 42*beta)*Power(mu,2))*tt + (-11 - 21*beta - 18*(11 + 3*beta)*Power(mu,2) + 5*(77 + 39*beta)*Power(mu,4))*Power(tt,2)))/(385.*Power(ss,2))

def S_12_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,2)*(-1 + 3*Power(mu,2)))/35.

def S_23_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-4 + 12*Power(mu,2) - 8*mu*(-3 + 5*Power(mu,2))*tt + (3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(tt,2)))/(35.*Power(ss,2))

def S_31_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,2)*(2 - 4*mu*tt + (-1 + 3*Power(mu,2))*Power(tt,2)))/(35.*Power(ss,2))

def T_12_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(4*mu*(143 + 78*beta - 45*Power(beta,2) + (715 + 3*beta*(364 + 155*beta))*Power(mu,2)) + (715 + 754*beta + 111*Power(beta,2) + 2*(429 + beta*(-182 + 69*beta))*Power(mu,2) - (5005 + 3*beta*(1690 + 643*beta))*Power(mu,4))*tt))/(10010.*tt) + (Power(beta,2)*(1144*mu*(-3*mu + (-2 + 5*Power(mu,2))*tt) + 3*Power(beta,2)*(37 + 46*Power(mu,2) - 643*Power(mu,4) + 5*mu*(-37 + 34*Power(mu,2) + 115*Power(mu,4))*tt) + 13*beta*(53 - 218*Power(mu,2) - 195*Power(mu,4) + mu*(-201 + 386*Power(mu,2) + 175*Power(mu,4))*tt)))/10010.

def T_23_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(572*mu*(-3 + 5*Power(mu,2))*Power(ss,4) + 26*Power(ss,2)*(8*mu*(11 + 3*beta*(-2 + 7*Power(mu,2))) + (88 - 264*Power(mu,2) - 6*beta*(-7 - 18*Power(mu,2) + 65*Power(mu,4)))*tt + mu*(44*(-3 + 5*Power(mu,2)) + 5*beta*(-21 + 10*Power(mu,2) + 35*Power(mu,4)))*Power(tt,2)) + beta*(60*mu*(26 + beta*(-3 + 31*Power(mu,2))) + (832 - 7072*Power(mu,2) + beta*(444 + 552*Power(mu,2) - 7716*Power(mu,4)))*tt + 18*mu*(104*(-2 + 7*Power(mu,2)) + 5*beta*(-37 + 34*Power(mu,2) + 115*Power(mu,4)))*Power(tt,2) - 12*(13*(-7 - 18*Power(mu,2) + 65*Power(mu,4)) + 5*beta*(-9 - 84*Power(mu,2) + 135*Power(mu,4) + 70*Power(mu,6)))*Power(tt,3) + 5*mu*(21*beta*(-9 - 10*Power(mu,2) + 35*Power(mu,4)) + 13*(-21 + 10*Power(mu,2) + 35*Power(mu,4)))*Power(tt,4))))/(10010.*Power(ss,4)*tt) + (Power(beta,2)*(-1144 + 416*beta + 111*Power(beta,2) - 3536*beta*Power(mu,2) + 138*Power(beta,2)*Power(mu,2) - 1929*Power(beta,2)*Power(mu,4) + 3*mu*(1144 + 624*beta*(-2 + 7*Power(mu,2)) + 15*Power(beta,2)*(-37 + 34*Power(mu,2) + 115*Power(mu,4)))*tt - 3*(572*(-1 + 3*Power(mu,2)) + 78*beta*(-7 - 18*Power(mu,2) + 65*Power(mu,4)) + 15*Power(beta,2)*(-9 - 84*Power(mu,2) + 135*Power(mu,4) + 70*Power(mu,6)))*Power(tt,2) + mu*(572*(-3 + 5*Power(mu,2)) + 105*Power(beta,2)*(-9 - 10*Power(mu,2) + 35*Power(mu,4)) + 130*beta*(-21 + 10*Power(mu,2) + 35*Power(mu,4)))*Power(tt,3) + 13*Power(ss,2)*(88*(1 - 3*Power(mu,2) + mu*(-3 + 5*Power(mu,2))*tt) + beta*(21 + 54*Power(mu,2) - 195*Power(mu,4) + 5*mu*(-21 + 10*Power(mu,2) + 35*Power(mu,4))*tt))))/(10010.*Power(ss,4))

def T_31_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*tt*(-8*mu*(143 + 390*beta + 210*Power(beta,2) + 13*(22 + 15*beta)*Power(ss,2)) + 4*(429*(-1 + 3*Power(mu,2)) + 156*beta*(-2 + 17*Power(mu,2)) + 45*Power(beta,2)*(-1 + 29*Power(mu,2)) + 26*(-11 + 33*Power(mu,2) + beta*(-2 + 17*Power(mu,2)))*Power(ss,2))*tt - 12*mu*(143*(-3 + 5*Power(mu,2)) + 156*beta*(-2 + 7*Power(mu,2)) + 15*Power(beta,2)*(-3 + 31*Power(mu,2)))*Power(tt,2) + (143*(3 - 30*Power(mu,2) + 35*Power(mu,4)) + 78*beta*(-7 - 18*Power(mu,2) + 65*Power(mu,4)) + 3*Power(beta,2)*(-37 - 46*Power(mu,2) + 643*Power(mu,4)))*Power(tt,3)))/(10010.*Power(ss,4)) - (Power(beta,2)*tt*(1144*mu*Power(ss,4) + 104*Power(ss,2)*(2*(11 + 15*beta)*mu + (22 - 66*Power(mu,2) + beta*(8 - 68*Power(mu,2)))*tt + mu*(-33 + 55*Power(mu,2) + 6*beta*(-2 + 7*Power(mu,2)))*Power(tt,2)) + beta*(120*(13 + 14*beta)*mu - 16*(-52 + 442*Power(mu,2) + 15*beta*(-1 + 29*Power(mu,2)))*tt + 72*mu*(26*(-2 + 7*Power(mu,2)) + 5*beta*(-3 + 31*Power(mu,2)))*Power(tt,2) - 12*(13*(-7 - 18*Power(mu,2) + 65*Power(mu,4)) + beta*(-37 - 46*Power(mu,2) + 643*Power(mu,4)))*Power(tt,3) + 5*mu*(13*(-21 + 10*Power(mu,2) + 35*Power(mu,4)) + 3*beta*(-37 + 34*Power(mu,2) + 115*Power(mu,4)))*Power(tt,4))))/(10010.*Power(ss,4))

#==End l=4,m=0




#==Starting l=4,m=1


def R_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Sqrt(1 - Power(mu,2))*(40*beta*mu + 44*mu*(1 + Power(ss,2)) - (11 + 3*beta + 3*(55 + 39*beta)*Power(mu,2))*tt + 2*mu*(-11 + 3*beta + (77 + 37*beta)*Power(mu,2))*Power(tt,2)))/(77.*Sqrt(5)*Power(ss,2))

def S_12_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,2)*mu*Sqrt(1 - Power(mu,2)))/(7.*Sqrt(5))

def S_23_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Sqrt(1 - Power(mu,2))*(-1 + 2*mu*tt)*(-4*mu + (-3 + 7*Power(mu,2))*tt))/(7.*Sqrt(5)*Power(ss,2))

def S_31_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,2)*Sqrt(1 - Power(mu,2))*tt*(-1 + mu*tt))/(7.*Sqrt(5)*Power(ss,2))

def T_12_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*Sqrt(1 - Power(mu,2))*(-715 - 598*beta - 75*Power(beta,2) - 3*(715 + beta*(1014 + 395*beta))*Power(mu,2) + 4*(143 + beta*(338 + 75*beta))*mu*tt + 4*(1001 + beta*(962 + 345*beta))*Power(mu,3)*tt))/(4004.*Sqrt(5)*tt) + (Power(beta,2)*Sqrt(1 - Power(mu,2))*(-2*(572 + beta*(559 + 75*beta))*mu - 2*beta*(481 + 345*beta)*Power(mu,3) - (143 + 3*beta*(52 + 5*beta))*tt + 3*(715 + 27*beta*(26 + 5*beta))*Power(mu,2)*tt + 10*beta*(91 + 66*beta)*Power(mu,4)*tt))/(2002.*Sqrt(5))

def T_23_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Sqrt(1 - Power(mu,2))*(-4*mu*(5*beta*(104 + 15*beta + 69*beta*Power(mu,2)) + 13*(44 + 3*beta + 37*beta*Power(mu,2))*Power(ss,2)) + 2*(858 + 117*beta*(1 + 39*Power(mu,2)) + 45*Power(beta,2)*(-1 + 27*Power(mu,2) + 44*Power(mu,4)) + 13*(33*(-1 + 5*Power(mu,2)) + 5*beta*(-3 + 9*Power(mu,2) + 14*Power(mu,4)))*Power(ss,2))*tt - 12*mu*(286 + 26*beta*(3 + 37*Power(mu,2)) + 15*Power(beta,2)*(-3 + 31*Power(mu,2) + 14*Power(mu,4)))*Power(tt,2) + (429*(-1 + 5*Power(mu,2)) + 260*beta*(-3 + 9*Power(mu,2) + 14*Power(mu,4)) + 105*Power(beta,2)*(-3 + 3*Power(mu,2) + 28*Power(mu,4)))*Power(tt,3)))/(4004.*Sqrt(5)*Power(ss,4)) - (Power(beta,2)*Sqrt(1 - Power(mu,2))*(429*(-1 + 5*Power(mu,2))*Power(ss,4) + 26*Power(ss,2)*(44 + 3*beta*(1 + 39*Power(mu,2)) - 8*mu*(22 + beta*(3 + 37*Power(mu,2)))*tt + (33*(-1 + 5*Power(mu,2)) + 10*beta*(-3 + 9*Power(mu,2) + 14*Power(mu,4)))*Power(tt,2)) + beta*(520 + 15*beta*(5 + 79*Power(mu,2)) - 80*mu*(52 + 3*beta*(5 + 23*Power(mu,2)))*tt + 18*(13 + 507*Power(mu,2) + 10*beta*(-1 + 27*Power(mu,2) + 44*Power(mu,4)))*Power(tt,2) - 16*mu*(39 + 481*Power(mu,2) + 15*beta*(-3 + 31*Power(mu,2) + 14*Power(mu,4)))*Power(tt,3) + 5*(26*(-3 + 9*Power(mu,2) + 14*Power(mu,4)) + 21*beta*(-3 + 3*Power(mu,2) + 28*Power(mu,4)))*Power(tt,4))))/(4004.*Sqrt(5)*Power(ss,4)*tt)

def T_31_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*Sqrt(1 - Power(mu,2))*tt*(-4*(143 + 5*beta*(52 + 21*beta)) + 24*(143 + 5*beta*(52 + 21*beta))*mu*tt - 9*(-143 + 26*beta + 25*Power(beta,2) + (715 + beta*(1014 + 395*beta))*Power(mu,2))*Power(tt,2) + 4*mu*(-429 + 78*beta + 75*Power(beta,2) + (1001 + beta*(962 + 345*beta))*Power(mu,2))*Power(tt,3) + 104*(11 + 5*beta)*Power(ss,2)*(-1 + 2*mu*tt)))/(4004.*Sqrt(5)*Power(ss,4)) - (Power(beta,2)*Sqrt(1 - Power(mu,2))*tt*(286*Power(ss,4) + 13*Power(ss,2)*(44 + 40*beta - 16*(11 + 10*beta)*mu*tt + 3*(-11 + beta + 55*Power(mu,2) + 39*beta*Power(mu,2))*Power(tt,2)) + beta*(260 + 210*beta - 80*(26 + 21*beta)*mu*tt + 9*(13 + 507*Power(mu,2) + 5*beta*(5 + 79*Power(mu,2)))*Power(tt,2) - 8*mu*(39 + 481*Power(mu,2) + 15*beta*(5 + 23*Power(mu,2)))*Power(tt,3) + 5*(13*(-3 + 9*Power(mu,2) + 14*Power(mu,4)) + 3*beta*(-1 + 27*Power(mu,2) + 44*Power(mu,4)))*Power(tt,4))))/(2002.*Sqrt(5)*Power(ss,4))

#==End l=4,m=1




#==Starting l=4,m=2


def R_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.4)*Power(beta,2)*(-1 + Power(mu,2))*(11 + 6*beta + 11*Power(ss,2) - 6*(11 + 6*beta)*mu*tt + 5*beta*Power(tt,2) + (77 + 31*beta)*Power(mu,2)*Power(tt,2)))/(77.*Power(ss,2))

def S_12_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.4)*Power(beta,2)*(-1 + Power(mu,2)))/7.

def S_23_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.4)*Power(beta,2)*(-1 + Power(mu,2))*(1 - 6*mu*tt + (-1 + 7*Power(mu,2))*Power(tt,2)))/(7.*Power(ss,2))

def S_31_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.4)*Power(beta,2)*(-1 + Power(mu,2))*Power(tt,2))/(7.*Power(ss,2))

def T_12_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-1 + Power(mu,2))*(-3*(143 + 3*beta*(52 + 15*beta))*mu + (143 + 208*beta + 33*Power(beta,2) + (1001 + beta*(806 + 237*beta))*Power(mu,2))*tt))/(1001.*Sqrt(10)*tt) - (Power(beta,2)*(-1 + Power(mu,2))*(286*(-1 + 3*mu*tt) + beta*(-221 - 33*beta - (403 + 237*beta)*Power(mu,2) + (793 + 165*beta)*mu*tt + 5*(91 + 57*beta)*Power(mu,3)*tt)))/(1001.*Sqrt(10))

def T_23_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-1 + Power(mu,2))*(429*mu*Power(ss,4) + 26*Power(ss,2)*(18*beta*mu - 2*(11 + beta*(5 + 31*Power(mu,2)))*tt + mu*(33 + 5*beta*(5 + 7*Power(mu,2)))*Power(tt,2)) + beta*(135*beta*mu - 12*(26 + beta*(11 + 79*Power(mu,2)))*tt + 18*mu*(78 + beta*(55 + 95*Power(mu,2)))*Power(tt,2) - 4*(65 + 403*Power(mu,2) + 15*beta*(2 + 29*Power(mu,2) + 14*Power(mu,4)))*Power(tt,3) + 5*mu*(65 + 91*Power(mu,2) + 21*beta*(2 + 7*Power(mu,2)))*Power(tt,4))))/(1001.*Sqrt(10)*Power(ss,4)*tt) - (Power(beta,2)*(-1 + Power(mu,2))*(-3*beta*(52 + beta*(11 + 79*Power(mu,2))) + 9*beta*mu*(156 + beta*(55 + 95*Power(mu,2)))*tt - 3*(143 + 26*beta*(5 + 31*Power(mu,2)) + 15*Power(beta,2)*(2 + 29*Power(mu,2) + 14*Power(mu,4)))*Power(tt,2) + mu*(429 + 105*Power(beta,2)*(2 + 7*Power(mu,2)) + 130*beta*(5 + 7*Power(mu,2)))*Power(tt,3) + 13*Power(ss,2)*(-22 + 66*mu*tt + beta*(-5 - 31*Power(mu,2) + 5*mu*(5 + 7*Power(mu,2))*tt))))/(1001.*Sqrt(10)*Power(ss,4))

def T_31_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-1 + Power(mu,2))*Power(tt,2)*(429 + 9*beta*(52 + 15*beta) + 26*(11 + 3*beta)*Power(ss,2) - 9*(143 + 3*beta*(52 + 15*beta))*mu*tt + (-143 + 130*beta + 33*Power(beta,2) + (1001 + beta*(806 + 237*beta))*Power(mu,2))*Power(tt,2)))/(1001.*Sqrt(10)*Power(ss,4)) + (Power(beta,2)*(-1 + Power(mu,2))*Power(tt,2)*(26*(11 + 6*beta)*Power(ss,2)*(-2 + 3*mu*tt) + beta*(-12*(26 + 15*beta) + 54*(26 + 15*beta)*mu*tt - 4*(65 + 33*beta + (403 + 237*beta)*Power(mu,2))*Power(tt,2) + 5*mu*(65 + 33*beta + (91 + 57*beta)*Power(mu,2))*Power(tt,3))))/(1001.*Sqrt(10)*Power(ss,4))

#==End l=4,m=2




#==Starting l=4,m=3


def R_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(11 + 3*beta)*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(11.*Sqrt(35)*Power(ss,2))

def S_12_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def S_23_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(Sqrt(35)*Power(ss,2))

def S_31_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def T_12_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(143 + 78*beta + 15*Power(beta,2))*Power(1 - Power(mu,2),1.5)*(-1 + 4*mu*tt))/(572.*Sqrt(35)*tt) + (Power(beta,2)*Power(1 - Power(mu,2),1.5)*(-6*beta*(13 + 5*beta)*mu + (143 + beta*(104 + 15*beta + 10*(13 + 6*beta)*Power(mu,2)))*tt))/(286.*Sqrt(35))

def T_23_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(1 - Power(mu,2),1.5)*(-60*Power(beta,2)*mu + 18*beta*(13 + 5*beta*(1 + 4*Power(mu,2)))*tt - 36*beta*mu*(26 + 5*beta*(3 + 2*Power(mu,2)))*Power(tt,2) + (143 + 5*beta*(52 + 21*beta + 4*(26 + 21*beta)*Power(mu,2)))*Power(tt,3) + 26*Power(ss,2)*(-6*beta*mu + (11 + 5*beta*(1 + 2*Power(mu,2)))*tt)))/(572.*Sqrt(35)*Power(ss,4)) - (Power(beta,2)*Power(1 - Power(mu,2),1.5)*(143*Power(ss,4) + 26*Power(ss,2)*(3*beta - 24*beta*mu*tt + (11 + 10*beta*(1 + 2*Power(mu,2)))*Power(tt,2)) + beta*(15*beta - 240*beta*mu*tt + 18*(13 + 10*beta*(1 + 4*Power(mu,2)))*Power(tt,2) - 48*mu*(13 + 5*beta*(3 + 2*Power(mu,2)))*Power(tt,3) + 5*(26 + 21*beta + 4*(13 + 21*beta)*Power(mu,2))*Power(tt,4))))/(572.*Sqrt(35)*Power(ss,4)*tt)

def T_31_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(143 + 78*beta + 15*Power(beta,2))*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(-3 + 4*mu*tt))/(572.*Sqrt(35)*Power(ss,4)) - (Power(beta,2)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(13*(11 + 3*beta)*Power(ss,2) + beta*(9*(13 + 5*beta) - 24*(13 + 5*beta)*mu*tt + 5*(13 + 3*beta + 26*Power(mu,2) + 12*beta*Power(mu,2))*Power(tt,2))))/(286.*Sqrt(35)*Power(ss,4))

#==End l=4,m=3




#==Starting l=4,m=4


def R_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(11 + beta)*Power(-1 + Power(mu,2),2)*Power(tt,2))/(11.*Sqrt(70)*Power(ss,2))

def S_12_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def S_23_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(-1 + Power(mu,2),2)*Power(tt,2))/(Sqrt(70)*Power(ss,2))

def S_31_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def T_12_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(143 + 26*beta + 3*Power(beta,2))*Power(-1 + Power(mu,2),2))/(286.*Sqrt(70)) + (Power(beta,3)*(13 + 3*beta)*Power(-1 + Power(mu,2),2)*(-1 + 5*mu*tt))/(286.*Sqrt(70))

def T_23_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*Power(-1 + Power(mu,2),2)*(-12*beta + 90*beta*mu*tt - 4*(13 + 15*beta + 30*beta*Power(mu,2))*Power(tt,2) + 5*(13 + 21*beta)*mu*Power(tt,3) + 26*Power(ss,2)*(-2 + 5*mu*tt)))/(286.*Sqrt(70)*Power(ss,4)) + (Power(beta,3)*Power(-1 + Power(mu,2),2)*(-3*beta + 45*beta*mu*tt - 3*(26 + 15*beta + 30*beta*Power(mu,2))*Power(tt,2) + 5*(26 + 21*beta)*mu*Power(tt,3) + 13*Power(ss,2)*(-1 + 5*mu*tt)))/(286.*Sqrt(70)*Power(ss,4))

def T_31_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(143 + 26*beta + 3*Power(beta,2))*Power(-1 + Power(mu,2),2)*Power(tt,4))/(286.*Sqrt(70)*Power(ss,4)) - (Power(beta,3)*(13 + 3*beta)*Power(-1 + Power(mu,2),2)*Power(tt,4)*(-4 + 5*mu*tt))/(286.*Sqrt(70)*Power(ss,4))

#==End l=4,m=4




#==Starting l=6,m=0


def R_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Power(beta,3)*(-4 + 12*Power(mu,2) - 8*mu*(-3 + 5*Power(mu,2))*tt + (3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(tt,2)))/(231.*Power(ss,2))

def S_12_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(75 + 6*beta + 30*(3 + 10*beta)*Power(mu,2) - 5*(105 + 106*beta)*Power(mu,4) - 15*(9 + 2*beta)*mu*tt - 10*(45 + 46*beta)*Power(mu,3)*tt + 21*(45 + 34*beta)*Power(mu,5)*tt))/3465. - (2*Power(beta,3)*(6*(20 + 13*beta)*mu - 10*(30 + 19*beta)*Power(mu,3) + (-3*(-5 + beta) - 30*(12 + 5*beta)*Power(mu,2) + 5*(105 + 53*beta)*Power(mu,4))*tt))/(3465.*tt)

def T_23_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*(8*mu*(30 + beta*(-39 + 95*Power(mu,2)) + 30*(-3 + 5*Power(mu,2))*Power(ss,2)) - 8*(60*(-1 + 3*Power(mu,2)) + beta*(-6 - 300*Power(mu,2) + 530*Power(mu,4)) + 15*(3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(ss,2))*tt + 12*mu*(60*(-3 + 5*Power(mu,2)) + beta*(-30 - 460*Power(mu,2) + 714*Power(mu,4)) + 5*(15 - 70*Power(mu,2) + 63*Power(mu,4))*Power(ss,2))*Power(tt,2) - 4*(30*(3 - 30*Power(mu,2) + 35*Power(mu,4)) + beta*(85 - 435*Power(mu,2) - 945*Power(mu,4) + 1743*Power(mu,6)))*Power(tt,3) + mu*(30*(15 - 70*Power(mu,2) + 63*Power(mu,4)) + 7*beta*(85 - 315*Power(mu,2) + 63*Power(mu,4) + 231*Power(mu,6)))*Power(tt,4)))/(6930.*Power(ss,4)*tt) + (Power(beta,3)*(4*(60 - 180*Power(mu,2) + beta*(3 + 150*Power(mu,2) - 265*Power(mu,4))) + 12*mu*(60*(-3 + 5*Power(mu,2)) + beta*(-15 - 230*Power(mu,2) + 357*Power(mu,4)))*tt - 3*(60*(3 - 30*Power(mu,2) + 35*Power(mu,4)) + beta*(85 - 435*Power(mu,2) - 945*Power(mu,4) + 1743*Power(mu,6)))*Power(tt,2) + mu*(60*(15 - 70*Power(mu,2) + 63*Power(mu,4)) + 7*beta*(85 - 315*Power(mu,2) + 63*Power(mu,4) + 231*Power(mu,6)))*Power(tt,3) + 30*Power(ss,2)*(-3 + 30*Power(mu,2) - 35*Power(mu,4) + mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*tt)))/(6930.*Power(ss,4))

def T_31_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,3)*tt*((3*(-15 + beta)*Power(-1 + Power(mu,2),2)*Pi*Power(tt,3))/2. - 2*mu*Pi*(-1 + mu*tt)*(15*Power(ss,2) + 2*(15 + 14*beta)*Power(-1 + mu*tt,2)) - 3*(-1 + Power(mu,2))*Pi*tt*(5*Power(ss,2) + (30 + 13*beta)*(-1 + mu*tt)*(-1 + 2*mu*tt))))/(3465.*Pi*Power(ss,4)) - (4*Power(beta,3)*tt*((-3*(-15 + 2*beta)*Power(-1 + Power(mu,2),2)*Pi*Power(tt,3)*(-4 + 5*mu*tt))/4. + 2*mu*Pi*(30*Power(ss,2)*Power(-1 + mu*tt,2) + (15 + 28*beta)*Power(-1 + mu*tt,4)) + 2*(-1 + Power(mu,2))*Pi*tt*(15*Power(ss,2)*(-2 + 3*mu*tt) + (15 + 13*beta)*Power(-1 + mu*tt,2)*(-2 + 5*mu*tt))))/(3465.*Pi*Power(ss,4))

#==End l=6,m=0




#==Starting l=6,m=1


def R_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Sqrt(0.09523809523809523)*Power(beta,3)*Sqrt(1 - Power(mu,2))*(-1 + 2*mu*tt)*(-4*mu + (-3 + 7*Power(mu,2))*tt))/(33.*Power(ss,2))

def S_12_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.09523809523809523)*Power(beta,3)*Sqrt(1 - Power(mu,2))*(5 + 3*beta - 15*(5 + 3*beta)*Power(mu,2) - 4*(10 + 3*beta)*mu*tt + 4*(35 + 17*beta)*Power(mu,3)*tt))/(165.*tt) + (Power(beta,3)*Sqrt(1 - Power(mu,2))*(-8*mu*(5 - 6*beta + (35 + 34*beta)*Power(mu,2)) + 5*(-7 - 3*beta - 2*(5 + 9*beta)*Power(mu,2) + 7*(15 + 11*beta)*Power(mu,4))*tt))/(330.*Sqrt(42))

def T_23_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Sqrt(1 - Power(mu,2))*(-8*mu*(20 - 6*beta + 34*beta*Power(mu,2) + 5*(-3 + 7*Power(mu,2))*Power(ss,2)) + 5*(-36 + 180*Power(mu,2) + 3*beta*(-3 - 18*Power(mu,2) + 77*Power(mu,4)) + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(ss,2))*tt - 6*mu*(40*(-3 + 7*Power(mu,2)) + 3*beta*(-15 - 10*Power(mu,2) + 81*Power(mu,4)))*Power(tt,2) + (50*(1 - 14*Power(mu,2) + 21*Power(mu,4)) + 7*beta*(5 - 60*Power(mu,2) + 45*Power(mu,4) + 66*Power(mu,6)))*Power(tt,3)))/(330.*Sqrt(42)*Power(ss,4)) - (Power(beta,3)*Sqrt(1 - Power(mu,2))*(40 + 4*beta*(-3 + 45*Power(mu,2)) - 64*mu*(5 + beta*(-3 + 17*Power(mu,2)))*tt + 30*(-6 + 30*Power(mu,2) + beta*(-3 - 18*Power(mu,2) + 77*Power(mu,4)))*Power(tt,2) - 8*mu*(20*(-3 + 7*Power(mu,2)) + 3*beta*(-15 - 10*Power(mu,2) + 81*Power(mu,4)))*Power(tt,3) + (25*(1 - 14*Power(mu,2) + 21*Power(mu,4)) + 7*beta*(5 - 60*Power(mu,2) + 45*Power(mu,4) + 66*Power(mu,6)))*Power(tt,4) + 10*Power(ss,2)*(-6 + 30*Power(mu,2) - 16*mu*(-3 + 7*Power(mu,2))*tt + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,2))))/(330.*Sqrt(42)*Power(ss,4)*tt)

def T_31_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.09523809523809523)*Power(beta,3)*Sqrt(1 - Power(mu,2))*tt*(-2*(10 + 7*beta) + 12*(10 + 7*beta)*mu*tt + 9*(5 + beta - 5*(5 + 3*beta)*Power(mu,2))*Power(tt,2) + 4*mu*(-3*(5 + beta) + (35 + 17*beta)*Power(mu,2))*Power(tt,3) + 10*Power(ss,2)*(-1 + 2*mu*tt)))/(165.*Power(ss,4)) - (Power(beta,3)*Sqrt(1 - Power(mu,2))*tt*(8*(5 + 7*beta) - 64*(5 + 7*beta)*mu*tt + 36*(-5 + 25*Power(mu,2) + beta*(-2 + 30*Power(mu,2)))*Power(tt,2) - 32*mu*(-15 - 6*beta + 35*Power(mu,2) + 34*beta*Power(mu,2))*Power(tt,3) + 5*(5 - 70*Power(mu,2) + 105*Power(mu,4) + beta*(-3 - 18*Power(mu,2) + 77*Power(mu,4)))*Power(tt,4) + 20*Power(ss,2)*(4 - 16*mu*tt + 3*(-1 + 5*Power(mu,2))*Power(tt,2))))/(330.*Sqrt(42)*Power(ss,4))

#==End l=6,m=1




#==Starting l=6,m=2


def R_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-8*Power(beta,3)*(-1 + Power(mu,2))*(1 - 6*mu*tt + (-1 + 7*Power(mu,2))*Power(tt,2)))/(33.*Sqrt(105)*Power(ss,2))

def S_12_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,3)*(-1 + Power(mu,2))*(-3*(2 + beta)*mu + (-1 + 2*(7 + 3*beta)*Power(mu,2))*tt))/(33.*Sqrt(105)*tt) - (4*Power(beta,3)*(-1 + Power(mu,2))*(-1 - (7 + 6*beta)*Power(mu,2) + (mu + 5*(3 + 2*beta)*Power(mu,3))*tt))/(33.*Sqrt(105))

def T_23_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(-1 + Power(mu,2))*(48*mu*(beta + 2*Power(ss,2)) - 64*(1 + 6*beta*Power(mu,2) + (-1 + 7*Power(mu,2))*Power(ss,2))*tt + 32*mu*(9 - 5*Power(ss,2) + 15*Power(mu,2)*(2*beta + Power(ss,2)))*Power(tt,2) - 4*(-16 - 15*beta + 2*(56 + 15*beta)*Power(mu,2) + 225*beta*Power(mu,4))*Power(tt,3) + mu*(-5*(16 + 21*beta) + 30*(8 + 7*beta)*Power(mu,2) + 231*beta*Power(mu,4))*Power(tt,4)))/(132.*Sqrt(105)*Power(ss,4)*tt) - (Power(beta,3)*(-1 + Power(mu,2))*(-32 - 96*beta*Power(mu,2) + 96*mu*(3 + 5*beta*Power(mu,2))*tt + (96 - 672*Power(mu,2) - 45*beta*(-1 + 2*Power(mu,2) + 15*Power(mu,4)))*Power(tt,2) + mu*(-5*(32 + 21*beta) + 30*(16 + 7*beta)*Power(mu,2) + 231*beta*Power(mu,4))*Power(tt,3) + 16*Power(ss,2)*(1 - 7*Power(mu,2) + 5*mu*(-1 + 3*Power(mu,2))*tt)))/(132.*Sqrt(105)*Power(ss,4))

def T_31_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,3)*(-1 + Power(mu,2))*Power(tt,2)*(3*(2 + beta) + Power(ss,2) - 9*(2 + beta)*mu*tt + 2*(-1 + (7 + 3*beta)*Power(mu,2))*Power(tt,2)))/(33.*Sqrt(105)*Power(ss,4)) + (4*Power(beta,3)*(-1 + Power(mu,2))*Power(tt,2)*(-4*(1 + beta) + 18*(1 + beta)*mu*tt - 4*(-1 + (7 + 6*beta)*Power(mu,2))*Power(tt,2) + 5*mu*(-1 + (3 + 2*beta)*Power(mu,2))*Power(tt,3) + Power(ss,2)*(-4 + 6*mu*tt)))/(33.*Sqrt(105)*Power(ss,4))

#==End l=6,m=2




#==Starting l=6,m=3


def R_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,3)*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(11.*Sqrt(105)*Power(ss,2))

def S_12_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Power(beta,3)*(3 + beta)*Power(1 - Power(mu,2),1.5)*(-1 + 4*mu*tt))/(33.*Sqrt(105)*tt) + (Power(beta,3)*Power(1 - Power(mu,2),1.5)*(-16*(3 + 2*beta)*mu + (9 + 5*beta)*(1 + 15*Power(mu,2))*tt))/(132.*Sqrt(105))

def T_23_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Power(1 - Power(mu,2),1.5)*(-16*mu*(2*beta + 3*Power(ss,2)) + 3*(24 + 5*beta + 75*beta*Power(mu,2) + 5*(-1 + 9*Power(mu,2))*Power(ss,2))*tt - 6*mu*(48 + 5*beta*(3 + 13*Power(mu,2)))*Power(tt,2) + (30*(-1 + 9*Power(mu,2)) + 7*beta*(-3 + 21*Power(mu,2) + 22*Power(mu,4)))*Power(tt,3)))/(132.*Sqrt(105)*Power(ss,4)) - (Power(beta,3)*Power(1 - Power(mu,2),1.5)*(8*beta - 128*beta*mu*tt + 6*(12 + beta*(5 + 75*Power(mu,2)))*Power(tt,2) - 8*mu*(24 + 5*beta*(3 + 13*Power(mu,2)))*Power(tt,3) + (15*(-1 + 9*Power(mu,2)) + 7*beta*(-3 + 21*Power(mu,2) + 22*Power(mu,4)))*Power(tt,4) + 6*Power(ss,2)*(4 - 32*mu*tt + 5*(-1 + 9*Power(mu,2))*Power(tt,2))))/(132.*Sqrt(105)*Power(ss,4)*tt)

def T_31_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Power(beta,3)*(3 + beta)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(-3 + 4*mu*tt))/(33.*Sqrt(105)*Power(ss,4)) - (Power(beta,3)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(72 + 48*beta + 24*Power(ss,2) - 64*(3 + 2*beta)*mu*tt + 5*(-3 + beta + 27*Power(mu,2) + 15*beta*Power(mu,2))*Power(tt,2)))/(132.*Sqrt(105)*Power(ss,4))

#==End l=6,m=3




#==Starting l=6,m=4


def R_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.2857142857142857)*Power(beta,3)*Power(-1 + Power(mu,2),2)*Power(tt,2))/(33.*Power(ss,2))

def S_12_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.2857142857142857)*Power(beta,3)*(5 + beta)*Power(-1 + Power(mu,2),2))/165. + (Power(beta,3)*(5 + 2*beta)*Power(-1 + Power(mu,2),2)*(-1 + 5*mu*tt))/(165.*Sqrt(14))

def T_23_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*Power(-1 + Power(mu,2),2)*(-16*beta + 120*beta*mu*tt - 4*(10 + 9*beta + 51*beta*Power(mu,2))*Power(tt,2) + mu*(50 + 63*beta + 77*beta*Power(mu,2))*Power(tt,3) + 20*Power(ss,2)*(-2 + 5*mu*tt)))/(330.*Sqrt(14)*Power(ss,4)) + (Power(beta,3)*Power(-1 + Power(mu,2),2)*(-4*beta + 60*beta*mu*tt - 3*(20 + 9*beta + 51*beta*Power(mu,2))*Power(tt,2) + mu*(100 + 63*beta + 77*beta*Power(mu,2))*Power(tt,3) + 10*Power(ss,2)*(-1 + 5*mu*tt)))/(330.*Sqrt(14)*Power(ss,4))

def T_31_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.2857142857142857)*Power(beta,3)*(5 + beta)*Power(-1 + Power(mu,2),2)*Power(tt,4))/(165.*Power(ss,4)) - (Power(beta,3)*(5 + 2*beta)*Power(-1 + Power(mu,2),2)*Power(tt,4)*(-4 + 5*mu*tt))/(165.*Sqrt(14)*Power(ss,4))

#==End l=6,m=4




#==Starting l=6,m=5


def R_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def S_12_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(5 + beta)*Power(1 - Power(mu,2),2.5)*tt)/(60.*Sqrt(77))

def T_23_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*Power(1 - Power(mu,2),2.5)*tt*(6*beta + 10*Power(ss,2) - 24*beta*mu*tt + (5 + 7*beta + 14*beta*Power(mu,2))*Power(tt,2)))/(60.*Sqrt(77)*Power(ss,4)) + (Power(beta,3)*Power(1 - Power(mu,2),2.5)*tt*(3*beta + 5*Power(ss,2) - 18*beta*mu*tt + (10 + 7*beta + 14*beta*Power(mu,2))*Power(tt,2)))/(60.*Sqrt(77)*Power(ss,4))

def T_31_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*(5 + beta)*Power(1 - Power(mu,2),2.5)*Power(tt,5))/(60.*Sqrt(77)*Power(ss,4))

#==End l=6,m=5




#==Starting l=6,m=6


def R_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def S_12_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def T_23_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,4)*Power(-1 + Power(mu,2),3)*Power(tt,2)*(-4 + 7*mu*tt))/(60.*Sqrt(231)*Power(ss,4)) - (Power(beta,4)*Power(-1 + Power(mu,2),3)*Power(tt,2)*(-3 + 7*mu*tt))/(60.*Sqrt(231)*Power(ss,4))

def T_31_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=6,m=6




#==Starting l=8,m=0


def R_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-8*Power(beta,4)*(4*mu*(3 - 5*Power(mu,2)) + (3 - 30*Power(mu,2) + 35*Power(mu,4))*tt))/(6435.*tt) + (8*Power(beta,4)*(-3 + 30*Power(mu,2) - 35*Power(mu,4) + mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*tt))/6435.

def T_23_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (32*Power(beta,4)*(2*Power(mu,4)*Pi*Power(-1 + mu*tt,3) + ((-1 + Power(mu,2))*Pi*(6 - 54*Power(mu,2) + 30*mu*(-3 + 11*Power(mu,2))*tt - 15*(1 - 20*Power(mu,2) + 43*Power(mu,4))*Power(tt,2) + 7*mu*(5 - 40*Power(mu,2) + 59*Power(mu,4))*Power(tt,3)))/8.))/(6435.*Pi*Power(ss,4)) - (32*Power(beta,4)*(2*Power(mu,3)*Pi*Power(-1 + mu*tt,4) + ((-1 + Power(mu,2))*Pi*(24*mu + tt*(24 - 216*Power(mu,2) + 60*mu*(-3 + 11*Power(mu,2))*tt - 20*(1 - 20*Power(mu,2) + 43*Power(mu,4))*Power(tt,2) + 7*mu*(5 - 40*Power(mu,2) + 59*Power(mu,4))*Power(tt,3))))/8.))/(6435.*Pi*Power(ss,4)*tt)

def T_31_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-32*Power(beta,4)*tt*(2*mu*Pi*Power(-1 + mu*tt,4) + (3*Power(-1 + Power(mu,2),2)*Pi*Power(tt,3)*(-4 + 5*mu*tt))/4. + 2*(-1 + Power(mu,2))*Pi*tt*Power(-1 + mu*tt,2)*(-2 + 5*mu*tt)))/(6435.*Pi*Power(ss,4)) - (32*Power(beta,4)*tt*(2*mu*Pi*Power(-1 + mu*tt,3) + (3*(-1 + Power(mu,2))*Pi*tt*(-2 + tt + 3*mu*tt)*(-2 + (-1 + 3*mu)*tt))/4.))/(6435.*Pi*Power(ss,4))

#==End l=8,m=0




#==Starting l=8,m=1


def R_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(2)*Power(beta,4)*Sqrt(1 - Power(mu,2))*(3 - 15*Power(mu,2) + 4*mu*(-3 + 7*Power(mu,2))*tt))/(2145.*tt) + (Sqrt(2)*Power(beta,4)*Sqrt(1 - Power(mu,2))*(8*mu*(3 - 7*Power(mu,2)) + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*tt))/2145.

def T_23_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,4)*Sqrt(1 - Power(mu,2))*(64*mu*(3 - 7*Power(mu,2)) + 120*(1 - 14*Power(mu,2) + 21*Power(mu,4))*tt - 144*mu*(5 - 30*Power(mu,2) + 33*Power(mu,4))*Power(tt,2) + 7*(-5 + 135*Power(mu,2) - 495*Power(mu,4) + 429*Power(mu,6))*Power(tt,3)))/(8580.*Sqrt(2)*Power(ss,4)) - (Power(beta,4)*Sqrt(1 - Power(mu,2))*(-48 + 240*Power(mu,2) + 256*mu*(3 - 7*Power(mu,2))*tt + 240*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,2) - 192*mu*(5 - 30*Power(mu,2) + 33*Power(mu,4))*Power(tt,3) + 7*(-5 + 135*Power(mu,2) - 495*Power(mu,4) + 429*Power(mu,6))*Power(tt,4)))/(8580.*Sqrt(2)*Power(ss,4)*tt)

def T_31_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(2)*Power(beta,4)*Sqrt(1 - Power(mu,2))*tt*(-4 + tt*(24*mu + tt*(9 - 45*Power(mu,2) + 4*mu*(-3 + 7*Power(mu,2))*tt))))/(2145.*Power(ss,4)) - (Sqrt(2)*Power(beta,4)*Sqrt(1 - Power(mu,2))*tt*(8 + tt*(-64*mu + tt*(-36 + 180*Power(mu,2) + 32*mu*(3 - 7*Power(mu,2))*tt + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,2)))))/(2145.*Power(ss,4))

#==End l=8,m=1




#==Starting l=8,m=2


def R_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-8*Power(beta,4)*(-1 + Power(mu,2))*(1 - 7*Power(mu,2) + 5*mu*(-1 + 3*Power(mu,2))*tt))/(429.*Sqrt(35)) + (8*Power(beta,4)*(-1 + Power(mu,2))*(-3*mu + (-1 + 7*Power(mu,2))*tt))/(429.*Sqrt(35)*tt)

def T_23_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,4)*(-1 + Power(mu,2))*(16 - 112*Power(mu,2) + 240*mu*(-1 + 3*Power(mu,2))*tt - 45*(1 - 18*Power(mu,2) + 33*Power(mu,4))*Power(tt,2) + 7*mu*(15 - 110*Power(mu,2) + 143*Power(mu,4))*Power(tt,3)))/(858.*Sqrt(35)*Power(ss,4)) + (Power(beta,4)*(-1 + Power(mu,2))*(48*mu + tt*(64 - 448*Power(mu,2) + 480*mu*(-1 + 3*Power(mu,2))*tt - 60*(1 - 18*Power(mu,2) + 33*Power(mu,4))*Power(tt,2) + 7*mu*(15 - 110*Power(mu,2) + 143*Power(mu,4))*Power(tt,3))))/(858.*Sqrt(35)*Power(ss,4)*tt)

def T_31_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (8*Power(beta,4)*(-1 + Power(mu,2))*Power(tt,2)*(3 - 9*mu*tt + (-1 + 7*Power(mu,2))*Power(tt,2)))/(429.*Sqrt(35)*Power(ss,4)) + (8*Power(beta,4)*(-1 + Power(mu,2))*Power(tt,2)*(-4 + 18*mu*tt + (4 - 28*Power(mu,2))*Power(tt,2) + 5*mu*(-1 + 3*Power(mu,2))*Power(tt,3)))/(429.*Sqrt(35)*Power(ss,4))

#==End l=8,m=2




#==Starting l=8,m=3


def R_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.0017316017316017316)*Power(beta,4)*Power(1 - Power(mu,2),1.5)*(-1 + 4*mu*tt))/(39.*tt) + (Power(beta,4)*Power(1 - Power(mu,2),1.5)*(-16*mu + 5*(-1 + 9*Power(mu,2))*tt))/(39.*Sqrt(2310))

def T_23_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,4)*Power(1 - Power(mu,2),1.5)*(16 - 256*mu*tt + 120*(-1 + 9*Power(mu,2))*Power(tt,2) - 160*mu*(-3 + 11*Power(mu,2))*Power(tt,3) + 7*(3 - 66*Power(mu,2) + 143*Power(mu,4))*Power(tt,4)))/(156.*Sqrt(2310)*Power(ss,4)*tt) + (Power(beta,4)*Power(1 - Power(mu,2),1.5)*(-64*mu + tt*(-60 + 540*Power(mu,2) + 120*mu*(3 - 11*Power(mu,2))*tt + 7*(3 - 66*Power(mu,2) + 143*Power(mu,4))*Power(tt,2))))/(156.*Sqrt(2310)*Power(ss,4))

def T_31_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.0017316017316017316)*Power(beta,4)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(-3 + 4*mu*tt))/(39.*Power(ss,4)) - (Power(beta,4)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(24 - 64*mu*tt + 5*(-1 + 9*Power(mu,2))*Power(tt,2)))/(39.*Sqrt(2310)*Power(ss,4))

#==End l=8,m=3




#==Starting l=8,m=4


def R_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2))/195. + (2*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*(-1 + 5*mu*tt))/195.

def T_23_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*(-8 + 60*mu*tt + (12 - 132*Power(mu,2))*Power(tt,2) + 7*mu*(-3 + 13*Power(mu,2))*Power(tt,3)))/(195.*Power(ss,4)) + (Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*(-2 + 30*mu*tt + (9 - 99*Power(mu,2))*Power(tt,2) + 7*mu*(-3 + 13*Power(mu,2))*Power(tt,3)))/(195.*Power(ss,4))

def T_31_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*Power(tt,4))/(195.*Power(ss,4)) - (2*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*Power(tt,4)*(-4 + 5*mu*tt))/(195.*Power(ss,4))

#==End l=8,m=4




#==Starting l=8,m=5


def R_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,4)*Power(1 - Power(mu,2),2.5)*tt)/(15.*Sqrt(2002))

def T_23_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,4)*Power(1 - Power(mu,2),2.5)*tt*(24 - 96*mu*tt + (-7 + 91*Power(mu,2))*Power(tt,2)))/(60.*Sqrt(2002)*Power(ss,4)) + (Power(beta,4)*Power(1 - Power(mu,2),2.5)*tt*(12 - 72*mu*tt + (-7 + 91*Power(mu,2))*Power(tt,2)))/(60.*Sqrt(2002)*Power(ss,4))

def T_31_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,4)*Power(1 - Power(mu,2),2.5)*Power(tt,5))/(15.*Sqrt(2002)*Power(ss,4))

#==End l=8,m=5




#==Starting l=8,m=6


def R_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def T_23_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,4)*Power(-1 + Power(mu,2),3)*Power(tt,2)*(-4 + 7*mu*tt))/(30.*Sqrt(429)*Power(ss,4)) - (Power(beta,4)*Power(-1 + Power(mu,2),3)*Power(tt,2)*(-3 + 7*mu*tt))/(30.*Sqrt(429)*Power(ss,4))

def T_31_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=8,m=6




#==Starting l=8,m=7


def R_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def T_23_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def T_31_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=8,m=7




#==Starting l=8,m=8


def R_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_12_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_23_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def S_31_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_12_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_23_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

def T_31_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0.

#==End l=8,m=8


