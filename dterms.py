import math
import kernels as kl
Power=pow
Sqrt=math.sqrt
Pi=math.pi


#==Starting l=0,m=0


def DR_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (7*(10 + beta + 2*beta*Power(mu,2) + beta*(1 + 2*Power(mu,2)))*Power(ss,2) - 2*(35 + 42*beta + 9*Power(beta,2))*mu*tt - 12*Power(beta,2)*Power(mu,3)*tt + (35 + beta*(28 + 3*beta))*(1 + Power(tt,2)) + 2*beta*(7 + 6*beta)*Power(mu,2)*(1 + Power(tt,2)) + beta*(-2*(42 + 18*beta)*mu*tt - 24*beta*Power(mu,3)*tt + (28 + 6*beta)*(1 + Power(tt,2)) + 12*beta*Power(mu,2)*(1 + Power(tt,2)) + 2*(7 + 6*beta)*Power(mu,2)*(1 + Power(tt,2))))/(105.*Power(ss,2))

def DS_12_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (10 + beta + 2*beta*Power(mu,2) + beta*(1 + 2*Power(mu,2)))/15.

def DS_23_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (5 + beta + 2*beta*Power(mu,2) + 5*Power(ss,2) - 2*(5 + 3*beta)*mu*tt + (5 + 3*beta)*Power(tt,2) + beta*(1 + 2*Power(mu,2) - 6*mu*tt + 3*Power(tt,2)))/(15.*Power(ss,2))

def DS_31_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (5 + 3*beta + 5*Power(ss,2) - 2*(5 + 3*beta)*mu*tt + (5 + beta + 2*beta*Power(mu,2))*Power(tt,2) + beta*(3 - 6*mu*tt + (1 + 2*Power(mu,2))*Power(tt,2)))/(15.*Power(ss,2))

def DT_12_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(mu*(105 + 378*beta + 20*Power(beta,3)*(3 + 4*Power(mu,2)) + 27*Power(beta,2)*(11 + 4*Power(mu,2))) - 2*(105 + 81*Power(beta,2)*(1 + 4*Power(mu,2)) + 42*beta*(5 + 4*Power(mu,2)) + 4*Power(beta,3)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*tt + mu*(105 + 378*beta + 20*Power(beta,3)*(3 + 4*Power(mu,2)) + 27*Power(beta,2)*(11 + 4*Power(mu,2)))*Power(tt,2)))/(630.*tt)

def DT_23_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(-21*(5 + 6*beta)*mu*Power(ss,4) - 3*Power(ss,2)*(6*beta*mu*(14 + beta*(9 + 6*Power(mu,2))) - (-35 + 28*beta*(1 + 2*Power(mu,2)) + 27*Power(beta,2)*(1 + 4*Power(mu,2)))*tt + 5*(-7 + 9*Power(beta,2))*mu*Power(tt,2)) + beta*(-5*beta*mu*(27 + 4*beta*(3 + 4*Power(mu,2))) + 6*(-21 + 9*beta*(1 + 4*Power(mu,2)) + 2*Power(beta,2)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*tt - 6*mu*(-63 + 10*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,2) + 2*(-63*(1 + 2*Power(mu,2)) - 27*beta*(1 + 4*Power(mu,2)) + 10*Power(beta,2)*(1 + 6*Power(mu,2)))*Power(tt,3) + 9*(14 + 15*beta)*mu*Power(tt,4))))/(630.*Power(ss,4)*tt)

def DT_31_00(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*tt*(21*(5 + 6*beta)*mu*Power(ss,4) + 3*Power(ss,2)*(5*(-7 + 9*Power(beta,2))*mu - (-35 + 28*beta*(1 + 2*Power(mu,2)) + 27*Power(beta,2)*(1 + 4*Power(mu,2)))*tt + 6*beta*mu*(14 + beta*(9 + 6*Power(mu,2)))*Power(tt,2)) + beta*(-9*(14 + 15*beta)*mu + (126*(1 + 2*Power(mu,2)) + 54*beta*(1 + 4*Power(mu,2)) - 20*Power(beta,2)*(1 + 6*Power(mu,2)))*tt + 6*mu*(-63 + 10*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,2) - 6*(-21 + 9*beta*(1 + 4*Power(mu,2)) + 2*Power(beta,2)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,3) + 5*beta*mu*(27 + 4*beta*(3 + 4*Power(mu,2)))*Power(tt,4))))/(630.*Power(ss,4))

#==End l=0,m=0




#==Starting l=2,m=0


def DR_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(13 + beta + (11 + 9*beta)*Power(mu,2) + beta*(1 + 9*Power(mu,2)) + (1 + 11*Power(mu,2))*Power(ss,2) - 2*mu*(3*beta + 3*(5 + beta) + 7*beta*Power(mu,2) + (9 + 7*beta)*Power(mu,2))*tt + (-5 + (29 + 12*beta)*Power(mu,2) + 8*beta*Power(mu,4))*Power(tt,2)))/(21.*Power(ss,2)) + (14 + beta*(13 + beta + (11 + 9*beta)*Power(mu,2)) + (7 + beta + (21 + 11*beta)*Power(mu,2))*Power(ss,2) - 2*mu*(14 + 3*beta*(5 + beta) + beta*(9 + 7*beta)*Power(mu,2))*tt + (-7 - 5*beta + (21 + beta*(29 + 6*beta))*Power(mu,2) + 4*Power(beta,2)*Power(mu,4))*Power(tt,2))/(21.*Power(ss,2))

def DS_12_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(1 + 11*Power(mu,2)))/21. + (7 + beta + 21*Power(mu,2) + 11*beta*Power(mu,2))/21.

def DS_23_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(1 + 11*Power(mu,2) - 2*mu*(3 + 9*Power(mu,2))*tt + 6*(-1 + 3*Power(mu,2))*Power(tt,2)))/(21.*Power(ss,2)) + (14 + beta + 11*beta*Power(mu,2) + 7*(-1 + 3*Power(mu,2))*Power(ss,2) - 2*mu*(14 + 3*beta + 9*beta*Power(mu,2))*tt + (7 + 6*beta)*(-1 + 3*Power(mu,2))*Power(tt,2))/(21.*Power(ss,2))

def DS_31_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(12 - 24*mu*tt + (1 + 11*Power(mu,2))*Power(tt,2)))/(21.*Power(ss,2)) + (2*(7 + 6*beta) + 14*Power(ss,2) - 4*(7 + 6*beta)*mu*tt + (-7 + beta + 21*Power(mu,2) + 11*beta*Power(mu,2))*Power(tt,2))/(21.*Power(ss,2))

def DT_12_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(2*mu*(231 + 297*beta*(3 + Power(mu,2)) + 99*Power(beta,2)*(8 + 7*Power(mu,2)) + 10*Power(beta,3)*(15 + 41*Power(mu,2))) - (231*(1 + 3*Power(mu,2)) + 132*beta*(5 + 31*Power(mu,2)) + 297*Power(beta,2)*(1 + 15*Power(mu,2) + 4*Power(mu,4)) + 16*Power(beta,3)*(3 + 69*Power(mu,2) + 68*Power(mu,4)))*tt + 2*mu*(231 + 594*beta*(1 + Power(mu,2)) + 297*Power(beta,2)*(1 + 4*Power(mu,2)) + 20*Power(beta,3)*(3 + 19*Power(mu,2) + 6*Power(mu,4)))*Power(tt,2)))/(1386.*tt)

def DT_23_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(33*mu*(7 + beta*(3 + 9*Power(mu,2)))*Power(ss,4) + 33*Power(ss,2)*(3*beta*mu*(8 + beta*(3 + 7*Power(mu,2))) - (-7 + 9*Power(beta,2)*Power(mu,2)*(3 + 2*Power(mu,2)) + beta*(2 + 22*Power(mu,2)))*tt + mu*(-7 + 15*Power(beta,2)*Power(mu,2))*Power(tt,2)) + beta*(5*beta*mu*(99 + beta*(30 + 82*Power(mu,2))) - 3*(-132 + 33*beta*(1 + 9*Power(mu,2)) + 4*Power(beta,2)*(3 + 69*Power(mu,2) + 68*Power(mu,4)))*tt + 12*mu*(-99 + 5*Power(beta,2)*(3 + 19*Power(mu,2) + 6*Power(mu,4)))*Power(tt,2) + (198*beta*Power(mu,2)*(3 + 2*Power(mu,2)) + 99*(1 + 11*Power(mu,2)) - 10*Power(beta,2)*(-1 + 21*Power(mu,2) + 36*Power(mu,4)))*Power(tt,3) - 99*(mu + (3 + 5*beta)*Power(mu,3))*Power(tt,4))))/(693.*Power(ss,4)*tt)

def DT_31_20(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*tt*(66*(7 + 12*beta)*mu*Power(ss,4) + 33*Power(ss,2)*(2*(-7 + 15*Power(beta,2))*mu - (7 - 21*Power(mu,2) + 9*Power(beta,2)*(1 + 9*Power(mu,2)) + beta*(4 + 44*Power(mu,2)))*tt + 6*beta*mu*(2 + 6*Power(mu,2) + beta*(3 + 7*Power(mu,2)))*Power(tt,2)) + 2*beta*(-99*(4 + 5*beta)*mu + (99*beta*(1 + 9*Power(mu,2)) + 99*(1 + 11*Power(mu,2)) - 10*Power(beta,2)*(5 + 51*Power(mu,2)))*tt + 3*mu*(-99*(1 + 3*Power(mu,2)) + 10*Power(beta,2)*(15 + 41*Power(mu,2)))*Power(tt,2) - 6*(33 - 99*Power(mu,2) + 33*beta*Power(mu,2)*(3 + 2*Power(mu,2)) + 2*Power(beta,2)*(3 + 69*Power(mu,2) + 68*Power(mu,4)))*Power(tt,3) + 5*beta*mu*(99*Power(mu,2) + 4*beta*(3 + 19*Power(mu,2) + 6*Power(mu,4)))*Power(tt,4))))/(1386.*Power(ss,4))

#==End l=2,m=0




#==Starting l=2,m=1


def DR_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.6666666666666666)*beta*Sqrt(1 - Power(mu,2))*(mu*(9 + 10*beta + 9*Power(ss,2)) - 3*(6 + beta + 2*(3 + 2*beta)*Power(mu,2) + beta*(1 + 4*Power(mu,2)))*tt + mu*(27 + 2*beta*(6 + 4*Power(mu,2)))*Power(tt,2)))/(21.*Power(ss,2)) + (Sqrt(0.6666666666666666)*Sqrt(1 - Power(mu,2))*(mu*(beta*(9 + 5*beta) + 3*(7 + 3*beta)*Power(ss,2)) - 3*(7 + beta*(6 + beta + 2*(3 + 2*beta)*Power(mu,2)))*tt + mu*(21 + beta*(27 + beta*(6 + 4*Power(mu,2))))*Power(tt,2)))/(21.*Power(ss,2))

def DS_12_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(6)*beta*mu*Sqrt(1 - Power(mu,2)))/7. + (Sqrt(0.6666666666666666)*(7 + 3*beta)*mu*Sqrt(1 - Power(mu,2)))/7.

def DS_23_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.6666666666666666)*beta*Sqrt(1 - Power(mu,2))*(3*mu - (3 + 6*Power(mu,2))*tt + 6*mu*Power(tt,2)))/(7.*Power(ss,2)) + (Sqrt(0.6666666666666666)*Sqrt(1 - Power(mu,2))*(mu*(3*beta + 7*Power(ss,2)) - (7 + beta*(3 + 6*Power(mu,2)))*tt + (7 + 6*beta)*mu*Power(tt,2)))/(7.*Power(ss,2))

def DS_31_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(6)*beta*Sqrt(1 - Power(mu,2))*tt*(-1 + mu*tt))/(7.*Power(ss,2)) + (Sqrt(0.6666666666666666)*(7 + 3*beta)*Sqrt(1 - Power(mu,2))*tt*(-1 + mu*tt))/(7.*Power(ss,2))

def DT_12_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*Sqrt(1 - Power(mu,2))*(231 + 198*beta*(3 + 2*Power(mu,2)) + 60*Power(beta,3)*(1 + 6*Power(mu,2)) + 33*Power(beta,2)*(11 + 24*Power(mu,2)) - 2*mu*(231 + 1188*beta + 80*Power(beta,3)*(3 + 4*Power(mu,2)) + 99*Power(beta,2)*(11 + 4*Power(mu,2)))*tt + (231 + 363*Power(beta,2)*(1 + 4*Power(mu,2)) + 198*beta*(3 + 4*Power(mu,2)) + 20*Power(beta,3)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,2)))/(462.*Sqrt(6)*tt)

def DT_23_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*Sqrt(1 - Power(mu,2))*(11*(7 + 6*beta*(1 + 2*Power(mu,2)))*Power(ss,4) + 11*Power(ss,2)*(6*beta*(2 + beta + 4*beta*Power(mu,2)) - 12*beta*mu*(2 + beta*(3 + 2*Power(mu,2)))*tt + (-7 + 5*Power(beta,2)*(1 + 4*Power(mu,2)))*Power(tt,2)) + beta*(5*beta*(11 + 4*beta*(1 + 6*Power(mu,2))) - 20*beta*mu*(11 + 4*beta*(3 + 4*Power(mu,2)))*tt + 2*(-99 + 10*Power(beta,2)*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,2) - 4*mu*(-99 - 22*beta*(3 + 2*Power(mu,2)) + 10*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,3) - 11*(6 + 12*Power(mu,2) + 5*beta*(1 + 4*Power(mu,2)))*Power(tt,4))))/(154.*Sqrt(6)*Power(ss,4)*tt)

def DT_31_21(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*Sqrt(1 - Power(mu,2))*tt*(33*(7 + 6*beta)*Power(ss,4) + 33*Power(ss,2)*(-7 + 5*Power(beta,2) - 2*(-7 + 12*beta + 15*Power(beta,2))*mu*tt + 6*beta*(2 + beta + 4*Power(mu,2) + 4*beta*Power(mu,2))*Power(tt,2)) + beta*(-33*(6 + 5*beta) + 4*(297 + 165*beta - 70*Power(beta,2))*mu*tt + 18*(-33*(1 + 2*Power(mu,2)) + 10*Power(beta,2)*(1 + 6*Power(mu,2)))*Power(tt,2) - 24*mu*(-33 + 11*beta*(3 + 2*Power(mu,2)) + 10*Power(beta,2)*(3 + 4*Power(mu,2)))*Power(tt,3) + 5*beta*(33*(1 + 4*Power(mu,2)) + 4*beta*(3 + 24*Power(mu,2) + 8*Power(mu,4)))*Power(tt,4))))/(462.*Sqrt(6)*Power(ss,4))

#==End l=2,m=1




#==Starting l=2,m=2


def DR_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2))*(3 + 2*beta + 3*Power(ss,2) - 6*beta*mu*tt - 6*(3 + beta)*mu*tt + (21 + 2*beta*(2 + 4*Power(mu,2)))*Power(tt,2)))/(21.*Sqrt(6)*Power(ss,2)) - ((-1 + Power(mu,2))*(beta*(3 + beta) + 3*(7 + beta)*Power(ss,2) - 6*beta*(3 + beta)*mu*tt + (21 + 21*beta + Power(beta,2)*(2 + 4*Power(mu,2)))*Power(tt,2)))/(21.*Sqrt(6)*Power(ss,2))

def DS_12_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2)))/(7.*Sqrt(6)) - ((7 + beta)*(-1 + Power(mu,2)))/(7.*Sqrt(6))

def DS_23_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2))*(1 - 6*mu*tt + 6*Power(tt,2)))/(7.*Sqrt(6)*Power(ss,2)) - ((-1 + Power(mu,2))*(beta + 7*Power(ss,2) - 6*beta*mu*tt + (7 + 6*beta)*Power(tt,2)))/(7.*Sqrt(6)*Power(ss,2))

def DS_31_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2))*Power(tt,2))/(7.*Sqrt(6)*Power(ss,2)) - ((7 + beta)*(-1 + Power(mu,2))*Power(tt,2))/(7.*Sqrt(6)*Power(ss,2))

def DT_12_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(beta*(-1 + Power(mu,2))*(6*beta*(33 + 33*beta + 10*Power(beta,2))*mu - 3*(77 + 220*beta + 16*Power(beta,3)*(1 + 4*Power(mu,2)) + 33*Power(beta,2)*(3 + 4*Power(mu,2)))*tt + 4*beta*mu*(99 + 132*beta + 10*Power(beta,2)*(3 + 2*Power(mu,2)))*Power(tt,2)))/(462.*Sqrt(6)*tt)

def DT_23_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-1 + Power(mu,2))*(10*Power(beta,2)*mu + 33*mu*Power(ss,4) - beta*(11 + 12*beta*(1 + 4*Power(mu,2)))*tt + 20*Power(beta,2)*mu*(3 + 2*Power(mu,2))*Power(tt,2) + (33 + 22*beta*(1 + 2*Power(mu,2)) - 10*Power(beta,2)*(1 + 4*Power(mu,2)))*Power(tt,3) - 11*(3 + 5*beta)*mu*Power(tt,4) + 11*Power(ss,2)*(3*beta*mu - (2 + 3*beta + 6*beta*Power(mu,2))*tt + 5*beta*mu*Power(tt,2))))/(77.*Sqrt(6)*Power(ss,4)*tt)

def DT_31_22(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (beta*(-1 + Power(mu,2))*Power(tt,2)*(2*beta*(99 + 33*beta - 10*Power(beta,2) + 9*(-33 + 10*Power(beta,2))*mu*tt - 6*(-33 + 11*beta*(1 + 2*Power(mu,2)) + 6*Power(beta,2)*(1 + 4*Power(mu,2)))*Power(tt,2) + 5*beta*mu*(33 + 4*beta*(3 + 2*Power(mu,2)))*Power(tt,3)) + 33*Power(ss,2)*(7 + 4*beta*(-1 + 3*mu*tt) + Power(beta,2)*(-3 + 6*mu*tt))))/(462.*Sqrt(6)*Power(ss,4))

#==End l=2,m=2




#==Starting l=4,m=0


def DR_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-16 + 136*Power(mu,2) + 8*mu*(12 - 42*Power(mu,2))*tt + (-21 - 54*Power(mu,2) + 195*Power(mu,4))*Power(tt,2)))/(385.*Power(ss,2)) + (2*beta*(44 - 16*beta + 4*(33 + 34*beta)*Power(mu,2) + 44*(-1 + 3*Power(mu,2))*Power(ss,2) + 8*mu*(11 + 12*beta - (55 + 42*beta)*Power(mu,2))*tt + (-11 - 21*beta - 18*(11 + 3*beta)*Power(mu,2) + 5*(77 + 39*beta)*Power(mu,4))*Power(tt,2)))/(385.*Power(ss,2))

def DS_12_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (8*beta*(-1 + 3*Power(mu,2)))/35.

def DS_23_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*beta*(-4 + 12*Power(mu,2) - 8*mu*(-3 + 5*Power(mu,2))*tt + (3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(tt,2)))/(35.*Power(ss,2))

def DS_31_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (8*beta*(2 - 4*mu*tt + (-1 + 3*Power(mu,2))*Power(tt,2)))/(35.*Power(ss,2))

def DT_12_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(8*mu*(143*(1 + 5*Power(mu,2)) + 117*beta*(1 + 14*Power(mu,2)) + 30*Power(beta,2)*(-3 + 31*Power(mu,2))) + (-286*(-5 + 18*Power(mu,2) + 35*Power(mu,4)) - 117*beta*(-37 + 82*Power(mu,2) + 195*Power(mu,4)) - 24*Power(beta,2)*(-37 - 46*Power(mu,2) + 643*Power(mu,4)))*tt + mu*(2288*(-2 + 5*Power(mu,2)) + 60*Power(beta,2)*(-37 + 34*Power(mu,2) + 115*Power(mu,4)) + 39*beta*(-201 + 386*Power(mu,2) + 175*Power(mu,4)))*Power(tt,2)))/(10010.*tt)

def DT_23_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(8*mu*(15*beta*(39 + beta*(-6 + 62*Power(mu,2))) + 26*(22 + 9*beta*(-2 + 7*Power(mu,2)))*Power(ss,2) + 143*(-3 + 5*Power(mu,2))*Power(ss,4)) + (-4*(-572 + 156*beta*(-2 + 17*Power(mu,2)) + 9*Power(beta,2)*(-37 - 46*Power(mu,2) + 643*Power(mu,4))) - 13*(176*(-1 + 3*Power(mu,2)) + 27*beta*(-7 - 18*Power(mu,2) + 65*Power(mu,4)))*Power(ss,2))*tt + 3*mu*(-2288 + 60*Power(beta,2)*(-37 + 34*Power(mu,2) + 115*Power(mu,4)) + 65*beta*(-21 + 10*Power(mu,2) + 35*Power(mu,4))*Power(ss,2))*Power(tt,2) - 6*(572*(1 - 3*Power(mu,2)) - 39*beta*(-7 - 18*Power(mu,2) + 65*Power(mu,4)) + 10*Power(beta,2)*(-9 - 84*Power(mu,2) + 135*Power(mu,4) + 70*Power(mu,6)))*Power(tt,3) - 13*mu*(88*(-3 + 5*Power(mu,2)) + 15*beta*(-21 + 10*Power(mu,2) + 35*Power(mu,4)))*Power(tt,4)))/(10010.*Power(ss,4)*tt)

def DT_31_40(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*tt*(104*mu*(-1 + Power(ss,2))*(22 + 45*beta + 22*Power(ss,2)) - 8*(3*(143 - 429*Power(mu,2) + beta*(52 - 442*Power(mu,2)) + 10*Power(beta,2)*(-1 + 29*Power(mu,2))) + 13*(-22 + 66*Power(mu,2) + 9*beta*(-2 + 17*Power(mu,2)))*Power(ss,2))*tt + 8*mu*(3*(429 - 715*Power(mu,2) + 30*Power(beta,2)*(-3 + 31*Power(mu,2))) + 26*(-33 + 55*Power(mu,2) + 9*beta*(-2 + 7*Power(mu,2)))*Power(ss,2))*Power(tt,2) - 2*(-143*(3 - 30*Power(mu,2) + 35*Power(mu,4)) + 117*beta*(-7 - 18*Power(mu,2) + 65*Power(mu,4)) + 18*Power(beta,2)*(-37 - 46*Power(mu,2) + 643*Power(mu,4)))*Power(tt,3) + 15*beta*mu*(13*(-21 + 10*Power(mu,2) + 35*Power(mu,4)) + 4*beta*(-37 + 34*Power(mu,2) + 115*Power(mu,4)))*Power(tt,4)))/(10010.*Power(ss,4))

#==End l=4,m=0




#==Starting l=4,m=1


def DR_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Sqrt(1 - Power(mu,2))*(40*mu - (3 + 117*Power(mu,2))*tt + 2*mu*(3 + 37*Power(mu,2))*Power(tt,2)))/(77.*Sqrt(5)*Power(ss,2)) + (2*beta*Sqrt(1 - Power(mu,2))*(40*beta*mu + 44*mu*(1 + Power(ss,2)) - (11 + 3*beta + 3*(55 + 39*beta)*Power(mu,2))*tt + 2*mu*(-11 + 3*beta + (77 + 37*beta)*Power(mu,2))*Power(tt,2)))/(77.*Sqrt(5)*Power(ss,2))

def DS_12_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (8*beta*mu*Sqrt(1 - Power(mu,2)))/(7.*Sqrt(5))

def DS_23_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*beta*Sqrt(1 - Power(mu,2))*(-1 + 2*mu*tt)*(-4*mu + (-3 + 7*Power(mu,2))*tt))/(7.*Sqrt(5)*Power(ss,2))

def DS_31_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (8*beta*Sqrt(1 - Power(mu,2))*tt*(-1 + mu*tt))/(7.*Sqrt(5)*Power(ss,2))

def DT_12_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Sqrt(1 - Power(mu,2))*(715*(1 + 3*Power(mu,2)) + 30*Power(beta,2)*(5 + 79*Power(mu,2)) + 39*beta*(23 + 117*Power(mu,2)) - 2*mu*(286*(5 + 7*Power(mu,2)) + 120*Power(beta,2)*(5 + 23*Power(mu,2)) + 117*beta*(23 + 37*Power(mu,2)))*tt + 2*(143*(-1 + 15*Power(mu,2)) + 39*beta*(-6 + 81*Power(mu,2) + 35*Power(mu,4)) + 30*Power(beta,2)*(-1 + 27*Power(mu,2) + 44*Power(mu,4)))*Power(tt,2)))/(2002.*Sqrt(5)*tt)

def DT_23_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Sqrt(1 - Power(mu,2))*(-429*(-1 + 5*Power(mu,2))*Power(ss,4) + 3*(-10*beta*(26 + beta*(5 + 79*Power(mu,2))) + 40*beta*mu*(26 + 3*beta*(5 + 23*Power(mu,2)))*tt - 4*(-143 + 15*Power(beta,2)*(-1 + 27*Power(mu,2) + 44*Power(mu,4)))*Power(tt,2) + 4*mu*(-286 - 13*beta*(3 + 37*Power(mu,2)) + 10*Power(beta,2)*(-3 + 31*Power(mu,2) + 14*Power(mu,4)))*Power(tt,3) + 13*(-11 + 55*Power(mu,2) + 5*beta*(-3 + 9*Power(mu,2) + 14*Power(mu,4)))*Power(tt,4)) - 13*Power(ss,2)*(88 - 176*mu*tt + 3*beta*(3 + 117*Power(mu,2) - 6*mu*(3 + 37*Power(mu,2))*tt + 5*(-3 + 9*Power(mu,2) + 14*Power(mu,4))*Power(tt,2)))))/(2002.*Sqrt(5)*Power(ss,4)*tt)

def DT_31_41(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*Sqrt(1 - Power(mu,2))*tt*(-52*(11 + 15*beta) + 572*Power(ss,4) - 24*(-143 - 130*beta + 70*Power(beta,2))*mu*tt + 9*(143 - 715*Power(mu,2) + 10*Power(beta,2)*(5 + 79*Power(mu,2)))*Power(tt,2) - 4*mu*(143*(3 - 7*Power(mu,2)) + 90*Power(beta,2)*(5 + 23*Power(mu,2)) + 39*beta*(3 + 37*Power(mu,2)))*Power(tt,3) + 15*beta*(13*(-3 + 9*Power(mu,2) + 14*Power(mu,4)) + 4*beta*(-1 + 27*Power(mu,2) + 44*Power(mu,4)))*Power(tt,4) + 13*Power(ss,2)*(60*beta - 8*(22 + 45*beta)*mu*tt + (-66 + 9*beta + 330*Power(mu,2) + 351*beta*Power(mu,2))*Power(tt,2))))/(2002.*Sqrt(5)*Power(ss,4))

#==End l=4,m=1




#==Starting l=4,m=2


def DR_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Sqrt(0.4)*Power(beta,2)*(-1 + Power(mu,2))*(6 - 36*mu*tt + 5*Power(tt,2) + 31*Power(mu,2)*Power(tt,2)))/(77.*Power(ss,2)) - (2*Sqrt(0.4)*beta*(-1 + Power(mu,2))*(11 + 6*beta + 11*Power(ss,2) - 6*(11 + 6*beta)*mu*tt + 5*beta*Power(tt,2) + (77 + 31*beta)*Power(mu,2)*Power(tt,2)))/(77.*Power(ss,2))

def DS_12_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.4)*beta*(-1 + Power(mu,2)))/7.

def DS_23_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.4)*beta*(-1 + Power(mu,2))*(1 - 6*mu*tt + (-1 + 7*Power(mu,2))*Power(tt,2)))/(7.*Power(ss,2))

def DS_31_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.4)*beta*(-1 + Power(mu,2))*Power(tt,2))/(7.*Power(ss,2))

def DT_12_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(-1 + Power(mu,2))*(6*(143 + 234*beta + 90*Power(beta,2))*mu - (286*(3 + 7*Power(mu,2)) + 117*beta*(11 + 31*Power(mu,2)) + 24*Power(beta,2)*(11 + 79*Power(mu,2)))*tt + 3*mu*(572 + 20*Power(beta,2)*(11 + 19*Power(mu,2)) + 13*beta*(61 + 35*Power(mu,2)))*Power(tt,2)))/(1001.*Sqrt(10)*tt)

def DT_23_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*(-1 + Power(mu,2))*(-540*Power(beta,2)*mu - 858*mu*Power(ss,4) + 36*beta*(13 + beta*(11 + 79*Power(mu,2)))*tt - 180*Power(beta,2)*mu*(11 + 19*Power(mu,2))*Power(tt,2) + 6*(-143 - 13*beta*(5 + 31*Power(mu,2)) + 10*Power(beta,2)*(2 + 29*Power(mu,2) + 14*Power(mu,4)))*Power(tt,3) + 39*mu*(22 + 5*beta*(5 + 7*Power(mu,2)))*Power(tt,4) - 13*Power(ss,2)*(108*beta*mu - (44 + 9*beta*(5 + 31*Power(mu,2)))*tt + 15*beta*mu*(5 + 7*Power(mu,2))*Power(tt,2))))/(1001.*Sqrt(10)*Power(ss,4)*tt)

def DT_31_42(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*(-1 + Power(mu,2))*Power(tt,2)*(858 + 468*beta - 180*Power(beta,2) + 18*(-143 + 90*Power(beta,2))*mu*tt - 2*(143*(1 - 7*Power(mu,2)) + 39*beta*(5 + 31*Power(mu,2)) + 18*Power(beta,2)*(11 + 79*Power(mu,2)))*Power(tt,2) + 15*beta*mu*(65 + 91*Power(mu,2) + beta*(44 + 76*Power(mu,2)))*Power(tt,3) + 26*Power(ss,2)*(-22 + 66*mu*tt + 27*beta*(-1 + 2*mu*tt))))/(1001.*Sqrt(10)*Power(ss,4))

#==End l=4,m=2




#==Starting l=4,m=3


def DR_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (3*Power(beta,2)*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(11.*Sqrt(35)*Power(ss,2)) + (2*beta*(11 + 3*beta)*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(11.*Sqrt(35)*Power(ss,2))

def DS_12_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*beta*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(Sqrt(35)*Power(ss,2))

def DS_31_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(1 - Power(mu,2),1.5)*(143 + 117*beta + 30*Power(beta,2) - 2*(286 + 351*beta + 120*Power(beta,2))*mu*tt + (286 + 60*Power(beta,2)*(1 + 4*Power(mu,2)) + 78*beta*(4 + 5*Power(mu,2)))*Power(tt,2)))/(286.*Sqrt(35)*tt)

def DT_23_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(1 - Power(mu,2),1.5)*(-30*Power(beta,2) - 143*Power(ss,4) + 360*Power(beta,2)*mu*tt - 180*Power(beta,2)*(1 + 4*Power(mu,2))*Power(tt,2) + 12*beta*mu*(-39 + 30*beta + 20*beta*Power(mu,2))*Power(tt,3) + 13*(11 + 15*beta + 30*beta*Power(mu,2))*Power(tt,4) - 39*beta*Power(ss,2)*(3 - 18*mu*tt + 5*(1 + 2*Power(mu,2))*Power(tt,2))))/(286.*Sqrt(35)*Power(ss,4)*tt)

def DT_31_43(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(-429 + 90*Power(beta,2) + 13*(22 + 9*beta)*Power(ss,2) - 4*(-143 + 117*beta + 90*Power(beta,2))*mu*tt + 15*beta*(13 + 4*beta + 26*Power(mu,2) + 16*beta*Power(mu,2))*Power(tt,2)))/(286.*Sqrt(35)*Power(ss,4))

#==End l=4,m=3




#==Starting l=4,m=4


def DR_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(-1 + Power(mu,2),2)*Power(tt,2))/(11.*Sqrt(70)*Power(ss,2)) + (Sqrt(0.05714285714285714)*beta*(11 + beta)*Power(-1 + Power(mu,2),2)*Power(tt,2))/(11.*Power(ss,2))

def DS_12_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.05714285714285714)*beta*Power(-1 + Power(mu,2),2)*Power(tt,2))/Power(ss,2)

def DS_31_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,2)*Power(-1 + Power(mu,2),2)*(-286 + 39*beta*(-3 + 5*mu*tt) + 12*Power(beta,2)*(-2 + 5*mu*tt)))/(286.*Sqrt(70))

def DT_23_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (3*Power(beta,3)*Power(-1 + Power(mu,2),2)*(12*beta - 60*beta*mu*tt + (-26 + 20*beta + 40*beta*Power(mu,2))*Power(tt,2) + 65*mu*Power(tt,3) + Power(ss,2)*(39 - 65*mu*tt)))/(286.*Sqrt(70)*Power(ss,4))

def DT_31_44(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,2)*Power(-1 + Power(mu,2),2)*Power(tt,4)*(286 + 12*Power(beta,2)*(-3 + 5*mu*tt) + 39*beta*(-2 + 5*mu*tt)))/(286.*Sqrt(70)*Power(ss,4))

#==End l=4,m=4




#==Starting l=6,m=0


def DR_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Power(beta,2)*(-4 + 12*Power(mu,2) - 8*mu*(-3 + 5*Power(mu,2))*tt + (3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(tt,2)))/(77.*Power(ss,2))

def DS_12_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(8*mu*(-90 - 78*beta + 225*Power(mu,2) + 190*beta*Power(mu,2)) + (beta*(48 + 2400*Power(mu,2) - 4240*Power(mu,4)) + 135*(1 + 18*Power(mu,2) - 35*Power(mu,4)))*tt + mu*(135*(-3 - 10*Power(mu,2) + 21*Power(mu,4)) + 8*beta*(-15 - 230*Power(mu,2) + 357*Power(mu,4)))*Power(tt,2)))/(3465.*tt)

def DT_23_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(-8*mu*(45 - 78*beta + 190*beta*Power(mu,2) + 45*(-3 + 5*Power(mu,2))*Power(ss,2)) + 3*(8*(-15 + 45*Power(mu,2) + beta*(-3 - 150*Power(mu,2) + 265*Power(mu,4))) + 45*(3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(ss,2))*tt + (24*beta*mu*(15 + 230*Power(mu,2) - 357*Power(mu,4)) - 45*mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*Power(ss,2))*Power(tt,2) + 2*(-45*(3 - 30*Power(mu,2) + 35*Power(mu,4)) + beta*(85 - 435*Power(mu,2) - 945*Power(mu,4) + 1743*Power(mu,6)))*Power(tt,3) + 45*mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*Power(tt,4)))/(3465.*Power(ss,4)*tt)

def DT_31_60(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*tt*(360*mu*(-1 + Power(ss,2)) - 4*(90 - 52*beta - 270*Power(mu,2) + 276*beta*Power(mu,2) + 135*(-1 + 3*Power(mu,2))*Power(ss,2))*tt + 24*mu*(2*beta*(-39 + 95*Power(mu,2)) + 15*(-3 + 5*Power(mu,2))*Power(ss,2))*Power(tt,2) - 6*(15*(3 - 30*Power(mu,2) + 35*Power(mu,4)) + 4*beta*(-3 - 150*Power(mu,2) + 265*Power(mu,4)))*Power(tt,3) + mu*(45*(15 - 70*Power(mu,2) + 63*Power(mu,4)) + 8*beta*(-15 - 230*Power(mu,2) + 357*Power(mu,4)))*Power(tt,4)))/(3465.*Power(ss,4))

#==End l=6,m=0




#==Starting l=6,m=1


def DR_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Sqrt(0.09523809523809523)*Power(beta,2)*Sqrt(1 - Power(mu,2))*(-1 + 2*mu*tt)*(-4*mu + (-3 + 7*Power(mu,2))*tt))/(11.*Power(ss,2))

def DS_12_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Sqrt(1 - Power(mu,2))*(12*(5 + 4*beta)*(-1 + 15*Power(mu,2)) - 8*mu*(-45 - 48*beta + 315*Power(mu,2) + 272*beta*Power(mu,2))*tt + 5*(-21 - 30*Power(mu,2) + 315*Power(mu,4) + 4*beta*(-3 - 18*Power(mu,2) + 77*Power(mu,4)))*Power(tt,2)))/(330.*Sqrt(42)*tt)

def DT_23_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Sqrt(1 - Power(mu,2))*(-8*(5 + beta*(-2 + 30*Power(mu,2))) + 32*mu*(5 + beta*(-6 + 34*Power(mu,2)))*tt - 20*beta*(-3 - 18*Power(mu,2) + 77*Power(mu,4))*Power(tt,2) + 8*mu*(30 - 70*Power(mu,2) + beta*(-15 - 10*Power(mu,2) + 81*Power(mu,4)))*Power(tt,3) + 25*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,4) - 5*Power(ss,2)*(-12 + 60*Power(mu,2) - 24*mu*(-3 + 7*Power(mu,2))*tt + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,2))))/(110.*Sqrt(42)*Power(ss,4)*tt)

def DT_31_61(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*Sqrt(1 - Power(mu,2))*tt*(-120 + 32*(15 - 14*beta)*mu*tt + 144*beta*(-1 + 15*Power(mu,2))*Power(tt,2) - 48*mu*(-15 - 12*beta + 35*Power(mu,2) + 68*beta*Power(mu,2))*Power(tt,3) + 5*(15*(1 - 14*Power(mu,2) + 21*Power(mu,4)) + 4*beta*(-3 - 18*Power(mu,2) + 77*Power(mu,4)))*Power(tt,4) + 60*Power(ss,2)*(2 - 12*mu*tt + 3*(-1 + 5*Power(mu,2))*Power(tt,2))))/(330.*Sqrt(42)*Power(ss,4))

#==End l=6,m=1




#==Starting l=6,m=2


def DR_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-8*Power(beta,2)*(-1 + Power(mu,2))*(1 - 6*mu*tt + (-1 + 7*Power(mu,2))*Power(tt,2)))/(11.*Sqrt(105)*Power(ss,2))

def DS_12_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-4*Power(beta,3)*mu*(-1 + Power(mu,2))*(6*(3 + 2*beta) - 3*(21 + 16*beta)*mu*tt + (3 + 5*(9 + 8*beta)*Power(mu,2))*Power(tt,2)))/(33.*Sqrt(105)*tt)

def DT_23_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(-1 + Power(mu,2))*(8*mu*(2*beta + 3*Power(ss,2)) - 4*(2 + 24*beta*Power(mu,2) + 3*(-1 + 7*Power(mu,2))*Power(ss,2))*tt + 20*(8*beta*Power(mu,3) + mu*(-1 + 3*Power(mu,2))*Power(ss,2))*Power(tt,2) + (-8 + 56*Power(mu,2) - 5*beta*(-1 + 2*Power(mu,2) + 15*Power(mu,4)))*Power(tt,3) + (20*mu - 60*Power(mu,3))*Power(tt,4)))/(11.*Sqrt(105)*Power(ss,4)*tt)

def DT_31_62(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,3)*(-1 + Power(mu,2))*Power(tt,2)*(6 - 4*beta + 36*beta*mu*tt - 6*(-1 + (7 + 12*beta)*Power(mu,2))*Power(tt,2) + 5*mu*(-3 + (9 + 8*beta)*Power(mu,2))*Power(tt,3) + 9*Power(ss,2)*(-1 + 2*mu*tt)))/(33.*Sqrt(105)*Power(ss,4))

#==End l=6,m=2




#==Starting l=6,m=3


def DR_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Sqrt(0.08571428571428572)*Power(beta,2)*Power(1 - Power(mu,2),1.5)*tt*(-1 + 2*mu*tt))/(11.*Power(ss,2))

def DS_12_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Power(1 - Power(mu,2),1.5)*(8*(9 + 4*beta) - 16*(27 + 16*beta)*mu*tt + (27 + 20*beta)*(1 + 15*Power(mu,2))*Power(tt,2)))/(132.*Sqrt(105)*tt)

def DT_23_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Power(1 - Power(mu,2),1.5)*(-32*beta + 384*beta*mu*tt - 60*(beta + 15*beta*Power(mu,2))*Power(tt,2) + 8*mu*(-36 + 15*beta + 65*beta*Power(mu,2))*Power(tt,3) + 45*(-1 + 9*Power(mu,2))*Power(tt,4) - 9*Power(ss,2)*(8 - 48*mu*tt + 5*(-1 + 9*Power(mu,2))*Power(tt,2))))/(132.*Sqrt(105)*Power(ss,4)*tt)

def DT_31_63(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(96*beta + 72*Power(ss,2) - 96*(3 + 4*beta)*mu*tt + 5*(-9 + 4*beta + 81*Power(mu,2) + 60*beta*Power(mu,2))*Power(tt,2)))/(132.*Sqrt(105)*Power(ss,4))

#==End l=6,m=3




#==Starting l=6,m=4


def DR_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Sqrt(0.2857142857142857)*Power(beta,2)*Power(-1 + Power(mu,2),2)*Power(tt,2))/(11.*Power(ss,2))

def DS_12_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Power(-1 + Power(mu,2),2)*(-45 + 75*mu*tt + 8*beta*(-2 + 5*mu*tt)))/(165.*Sqrt(14))

def DT_23_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Power(-1 + Power(mu,2),2)*(8*beta - 40*beta*mu*tt + 2*(-5 + 3*beta + 17*beta*Power(mu,2))*Power(tt,2) + 25*mu*Power(tt,3) - 5*Power(ss,2)*(-3 + 5*mu*tt)))/(55.*Sqrt(14)*Power(ss,4))

def DT_31_64(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*Power(-1 + Power(mu,2),2)*Power(tt,4)*(-30 + 75*mu*tt + 8*beta*(-3 + 5*mu*tt)))/(165.*Sqrt(14)*Power(ss,4))

#==End l=6,m=4




#==Starting l=6,m=5


def DR_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*(15 + 4*beta)*Power(1 - Power(mu,2),2.5)*tt)/(60.*Sqrt(77))

def DT_23_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (Power(beta,3)*Power(1 - Power(mu,2),2.5)*tt*(-4*beta - 5*Power(ss,2) + 8*beta*mu*tt + 5*Power(tt,2)))/(20.*Sqrt(77)*Power(ss,4))

def DT_31_65(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,3)*(15 + 4*beta)*Power(1 - Power(mu,2),2.5)*Power(tt,5))/(60.*Sqrt(77)*Power(ss,4))

#==End l=6,m=5




#==Starting l=6,m=6


def DR_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_23_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return -(Power(beta,4)*Power(-1 + Power(mu,2),3)*Power(tt,2))/(15.*Sqrt(231)*Power(ss,4))

def DT_31_66(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=6,m=6




#==Starting l=8,m=0


def DR_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (32*Power(beta,4)*(4*mu*(-3 + 5*Power(mu,2)) + (-6 + 60*Power(mu,2) - 70*Power(mu,4))*tt + mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*Power(tt,2)))/(6435.*tt)

def DT_23_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (16*Power(beta,4)*(8*mu*(3 - 5*Power(mu,2)) + 6*(3 - 30*Power(mu,2) + 35*Power(mu,4))*tt - 6*mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*Power(tt,2) + (-5 + 105*Power(mu,2) - 315*Power(mu,4) + 231*Power(mu,6))*Power(tt,3)))/(6435.*Power(ss,4)*tt)

def DT_31_80(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-32*Power(beta,4)*Power(tt,2)*(4 - 12*Power(mu,2) + 12*mu*(-3 + 5*Power(mu,2))*tt - 3*(3 - 30*Power(mu,2) + 35*Power(mu,4))*Power(tt,2) + mu*(15 - 70*Power(mu,2) + 63*Power(mu,4))*Power(tt,3)))/(6435.*Power(ss,4))

#==End l=8,m=0




#==Starting l=8,m=1


def DR_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,4)*Sqrt(2 - 2*Power(mu,2))*(-6 + 30*Power(mu,2) - 16*mu*(-3 + 7*Power(mu,2))*tt + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,2)))/(2145.*tt)

def DT_23_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Power(beta,4)*Sqrt(2 - 2*Power(mu,2))*(2 - 10*Power(mu,2) + 8*mu*(-3 + 7*Power(mu,2))*tt - 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,2) + 2*mu*(5 - 30*Power(mu,2) + 33*Power(mu,4))*Power(tt,3)))/(715.*Power(ss,4)*tt)

def DT_31_81(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-4*Power(beta,4)*Sqrt(2 - 2*Power(mu,2))*Power(tt,2)*(-16*mu + 18*(-1 + 5*Power(mu,2))*tt - 24*mu*(-3 + 7*Power(mu,2))*Power(tt,2) + 5*(1 - 14*Power(mu,2) + 21*Power(mu,4))*Power(tt,3)))/(2145.*Power(ss,4))

#==End l=8,m=1




#==Starting l=8,m=2


def DR_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-32*Power(beta,4)*(-1 + Power(mu,2))*(3*mu + (2 - 14*Power(mu,2))*tt + 5*mu*(-1 + 3*Power(mu,2))*Power(tt,2)))/(429.*Sqrt(35)*tt)

def DT_23_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Power(beta,4)*(-1 + Power(mu,2))*(-16*mu + 16*(-1 + 7*Power(mu,2))*tt - 80*mu*(-1 + 3*Power(mu,2))*Power(tt,2) + 5*(1 - 18*Power(mu,2) + 33*Power(mu,4))*Power(tt,3)))/(143.*Sqrt(35)*Power(ss,4)*tt)

def DT_31_82(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (32*Power(beta,4)*(-1 + Power(mu,2))*Power(tt,2)*(-1 + 9*mu*tt + (3 - 21*Power(mu,2))*Power(tt,2) + 5*mu*(-1 + 3*Power(mu,2))*Power(tt,3)))/(429.*Sqrt(35)*Power(ss,4))

#==End l=8,m=2




#==Starting l=8,m=3


def DR_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Sqrt(0.0017316017316017316)*Power(beta,4)*Power(1 - Power(mu,2),1.5)*(4 - 32*mu*tt + 5*(-1 + 9*Power(mu,2))*Power(tt,2)))/(39.*tt)

def DT_23_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Sqrt(0.0017316017316017316)*Power(beta,4)*Power(1 - Power(mu,2),1.5)*(-4 + 48*mu*tt - 15*(-1 + 9*Power(mu,2))*Power(tt,2) + 10*mu*(-3 + 11*Power(mu,2))*Power(tt,3)))/(39.*Power(ss,4)*tt)

def DT_31_83(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.0017316017316017316)*Power(beta,4)*Power(1 - Power(mu,2),1.5)*Power(tt,3)*(12 - 48*mu*tt + 5*(-1 + 9*Power(mu,2))*Power(tt,2)))/(39.*Power(ss,4))

#==End l=8,m=3




#==Starting l=8,m=4


def DR_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (8*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*(-2 + 5*mu*tt))/195.

def DT_23_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (4*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*(2 - 10*mu*tt + (-1 + 11*Power(mu,2))*Power(tt,2)))/(65.*Power(ss,4))

def DT_31_84(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-8*Sqrt(0.025974025974025976)*Power(beta,4)*Power(-1 + Power(mu,2),2)*Power(tt,4)*(-3 + 5*mu*tt))/(195.*Power(ss,4))

#==End l=8,m=4




#==Starting l=8,m=5


def DR_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Sqrt(0.001998001998001998)*Power(beta,4)*Power(1 - Power(mu,2),2.5)*tt)/15.

def DT_23_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (2*Sqrt(0.001998001998001998)*Power(beta,4)*Power(1 - Power(mu,2),2.5)*tt*(-1 + 2*mu*tt))/(5.*Power(ss,4))

def DT_31_85(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Sqrt(0.001998001998001998)*Power(beta,4)*Power(1 - Power(mu,2),2.5)*Power(tt,5))/(15.*Power(ss,4))

#==End l=8,m=5




#==Starting l=8,m=6


def DR_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_23_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return (-2*Power(beta,4)*Power(-1 + Power(mu,2),3)*Power(tt,2))/(15.*Sqrt(429)*Power(ss,4))

def DT_31_86(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=8,m=6




#==Starting l=8,m=7


def DR_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_23_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_31_87(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=8,m=7




#==Starting l=8,m=8


def DR_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_12_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_23_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DS_31_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_12_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_23_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

def DT_31_88(mu,tt,beta):
	ss=kl.sfac(mu,tt)
	return 0

#==End l=8,m=8


