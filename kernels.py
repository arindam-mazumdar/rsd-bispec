import math

"Defeinition of s in terms of mu and t"
#=======================================
def sfac(mu,tt):
	return math.sqrt(1-2*mu*tt+pow(tt,2.))

"Definition of G interms of mu and t"
#=======================================
def G12(mu,tt):
	return (6 - (7*mu)/tt - 7*tt*mu + 8*pow(mu,2))/14.

def G23(mu,tt):
	ss=sfac(mu,tt)#math.sqrt(1-2*mu*tt+pow(tt,2))
	return -(tt - 7*mu + 6*tt*pow(mu,2))/(14.*pow(ss,2)*tt)

def G31(mu,tt):
	ss=sfac(mu,tt)# math.sqrt(1-2*mu*tt+pow(tt,2))
	return (pow(tt,2)*(-1 + 7*tt*mu - 6*pow(mu,2)))/(14.*pow(ss,2))
"Definition of F in terms of mu and t"
#=======================================
def F12(mu,tt):
	return (10 - (7*mu)/tt - 7*tt*mu + 4*pow(mu,2))/14.

def F23(mu,tt):
        ss=sfac(mu,tt)# math.sqrt(1-2*mu*tt+pow(tt,2))
        return (7*mu + tt*(3 - 10*pow(mu,2)))/(14.*tt*ss*ss)

def F31(mu,tt):
        ss=sfac(mu,tt)# math.sqrt(1-2*mu*tt+pow(tt,2))
        return (pow(tt,2)*(3 + 7*tt*mu - 10*pow(mu,2)))/(14.*ss*ss)


"Definition of Delta-G in terms of mu and t"
#=======================================
def dG12(mu,tt):
	return (-2*(-1 + pow(mu,2)))/7.

def dG23(mu,tt):
	ss=sfac(mu,tt)
	return (-2*(-1 + pow(mu,2)))/(7.*pow(ss,2))

def dG31(mu,tt):
	ss=sfac(mu,tt)
	return (-2*pow(tt,2)*(-1 + pow(mu,2)))/(7.*pow(ss,2))
