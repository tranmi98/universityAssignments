import scipy.io 
import matplotlib as mpl
from scipy.optimize import fsolve
# Import math Library
import math 
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, solve_ivp

#######################################################################
#=====================================================================
# Nonlinear Dynamical Systems Homework 2
# Minh-Phuong Tran	51091600
# 
#
# If you run this file, it will produce the figures for question 3. 
# The other figures can be obtained by uncommenting the corresponding
# lines in the script area at the bottom of this file.
#=====================================================================
######################################################################

########################################################
# Functions in the ode
########################################################
def f(x,*args):
	r = args[1] 
	k = args[2]

	return r*x*(1-x/k)-x*x/(1+x*x)
	
def f2(x,*args):
	r = args[1]
	k = args[2]

	return r*(1-x/k)-x/(1+x*x)

def g(x,*args):
	r = args[1] #or args[0] if use rootsdiag
	k = args[2]

	return -r*x*(1-x/k)+x*x/(1+x*x)

def f_rootsdiag(x,*args):
	r = args[0] #or args[0] if use rootsdiag
	k = args[1]

	return r*x*(1-x/k)-x*x/(1+x*x)

def f2_rootsdiag(x,*args):
	r = args[0]
	k = args[1]

	return r*(1-x/k)-x/(1+x*x)



#################################################
# Plot intersections between \frac{x}{1+ x^2} and r(1 - \frac{x}{k})
##################################################
def plotIntersect(r,k,x):

	y1 = r*(1-x/k)
	y2 = x/(1+x*x)
	ax = plt.gca()
	ax.plot(x,y1, label = r'$f_1(x) = r(1 - \frac{x}{k})$')
	ax.plot(x,y2, label = r'$f_2(x) = \frac{x}{1+ x^2}$')

	idx = np.argwhere(np.diff(np.sign(y2 - y1))).flatten()
	plt.plot(x[idx], y1[idx], 'ro')
	txt = ["A","B","C"]
	for i in range(len(idx)):
	   	ax.annotate(txt[i], (x[idx[i]], y1[idx[i]]))

	plt.ylim(bottom = 0)
	plt.xlim(left = 0)
	plt.xlabel('x')
	ax.legend(loc='upper right')
	plt.show()


###################################################
# Plot bifurcation diagram as (k(x), r(x)) with 3 areas
###################################################
def plotbifurcurves(x):
	r_x = 2*np.power(x,3)/ np.power((1+x*x),2)
	k_x = 2*np.power(x,3)/ (x*x-1)
	plt.plot(k_x,r_x)
	plt.xlim(right = 80)
	plt.xlabel('k')
	plt.ylabel('r')
	plt.text(40, 0.3, '3 intersections')
	plt.text(5, 0.1, '1 intersection (i)')
	plt.text(40, 0.6, '1 intersection (ii)')
	plt.title('Bifurcation curves')
	plt.show()


###############################################################
# Roots of f or f2 in function of k
###############################################################
def rootsdiag(r):
	k = np.linspace(2,16,num = 200)
	x0 = [0.1,2,5,20]
	for i in range(len(k)):
		for j in range(4):
			(x,info,ier,msg) = fsolve(f_rootsdiag,x0[j], args = (r,k[i]),full_output=True)
			if ier == 1:
				plt.scatter(k[i],x)


	plt.xlabel('k')
	plt.ylabel('x')
	plt.title(r'Roots of f(x) = ($0.58x(1-\frac{x}{k}) - \frac{x^2}{1+x^2}$) in function of k')
	plt.show()

################################################################
# Different trajectories for different initial values
################################################################
def solvediffeq(t,r,k):
	y0 = np.linspace(0,5,30)
	for i in range(len(y0)):
		sol = odeint(f,y0[i],t,args = (r,k))
		plt.plot(t,sol)

	plt.xlabel('t')
	plt.ylabel('x')
	plt.title(r'Trajectories for r = 0.5 and k = 15 with initial value $x_0 \in [0,5]$')
	plt.show()

#################################################################
# Question 3 :
# Evolution of equilibrium and hysteresis 
# Question 4 :
# Stable part of bifurcation diagram
#################################################################
def solvevaryk(f,t,tau,r,k):
	y0 = 1.1
	ends = np.zeros((len(k),1))
	for i in range(len(k)):
		sol = odeint(f,y0,t,args = (r,k[i]))
		ends[i] = sol[-1]
		y0 = sol[-1]

	plt.scatter(k,ends)
	plt.show()
	plt.quiver(k[:-1], ends[:-1], k[1:]-k[:-1], ends[1:]-ends[:-1], scale_units='xy', angles='xy', scale=1)
	plt.xlabel('k')
	plt.ylabel('x')
	plt.title('Evolution of x at equilibrium with varying sinusoidal k')
	plt.show()

###################################################################
# Question 3 : 
# Trajectory of x(t) when k varies like a sinusoidal
####################################################################
def trajsolvevaryk(f,t,tau,r):
	k = np.sin(tau)*10+5
	y0 = 1.1
	t_step = len(t)//len(k)
	for i in range(len(k)-1):
		sol = odeint(f,y0,np.linspace(0,t[t_step],t_step),args = (r,k[i]))
		y0 = sol[-1]
		plt.plot(t[i*t_step:(i+1)*t_step],sol)

	plt.xlabel('t')
	plt.ylabel('x')
	plt.title('Evolution of x(t) with varying sinusoidal k')
	plt.show()

###################################################################
# Question 5 :
# Unstable part of bifurcation diagram
##################################################################
def solvevaryk2(t,tau,r):
	k = np.sin(tau)+6.15
	#plt.plot(t,k)
	#plt.show()
	y0 = 1.1
	ends = np.zeros((len(k),1))
	for i in range(len(k)):
		sol = odeint(g,y0,t,args = (r,k[i]))
		ends[i] = sol[-1]
		y0 = sol[-1]

	plt.scatter(k,ends)
	plt.xlabel('k')
	plt.ylabel('x')
	plt.title('Unstable part of the bifurcation diagram')
	plt.show()
####################################################################
# Question 5 :
# Both parts of bifurcation diagram
#####################################################################
def solvevaryk3(t,tau,r):
	k2 = 1.7*np.sin(tau)+6.15
	y0 = 1.1
	ends = np.zeros((len(k2),1))
	for i in range(len(k2)):
		sol = odeint(g,y0,t,args = (r,k2[i]))
		ends[i] = sol[-1]
		y0 = sol[-1]

	k = np.sin(tau)*10+5
	plt.scatter(k2,ends)
	y0 = 1.1
	ends = np.zeros((len(k),1))
	for i in range(len(k)):
		sol = odeint(f,y0,t,args = (r,k[i]))
		ends[i] = sol[-1]
		y0 = sol[-1]

	plt.scatter(k,ends)
	plt.xlabel('k')
	plt.ylabel('x')
	plt.title('Both parts of the bifurcation diagram')
	plt.show()


r = 0.58

#==============================================================
#SCRIPTS
#==============================================================
#--------------------------------------------------------------
# Plot 3 intersections
#--------------------------------------------------------------
#r = 0.4
#k = 20
#x = np.linspace(0,2*k,num = 2000)
#plotIntersect(r,k,x)

#-------------------------------------------------------------
#Plot roots of polynomial f for varying k
#--------------------------------------------------------------
#r = 0.58
#rootsdiag(r)

#--------------------------------------------------------------
# Plot Bifurcation diagram (k(x), r(x))
#--------------------------------------------------------------
#x = np.linspace(1.01,40,num = 2000)
#plotbifurcurves(x)

#-----------------------------------------------------------------
# Plot different trajectories for different initial conditions for
# r = 0.58
# k = 7
#------------------------------------------------------------------
#r= 0.58
#k = 7
#t = np.linspace(0,80,1000)
#solvediffeq(t,r,k)


#-----------------------------------------------------------------
# QUESTION 3 & 4
# Plot evolution of equilibrium with varying k
# r = 0.58
# x_0 = 1.1
#-----------------------------------------------------------------
r = 0.58
tau = np.linspace(0,math.pi,num = 50)
t = np.linspace(0,50,200)
k = np.sin(tau)*10+5
#plt.plot(tau,k)
#plt.xlabel('t')
#plt.ylabel('k')
#plt.title('Sinusoidal k(t) = 10sin(t)+5')
#plt.show()
solvevaryk(f,t,tau,r,k)


#-----------------------------------------------------------------
# QUESTION 3
# Plot trajectory of x(t) with varying k
# r = 0.58
# x_0 = 1.1
#-----------------------------------------------------------------
r = 0.58
tau = np.linspace(0,math.pi,num = 50)
t = np.linspace(0,9000,9000)
trajsolvevaryk(f,t,tau,r)

#-----------------------------------------------------------------
# QUESTION 5
# Unstable part of bifurcation diagram
# r = 0.58
# x_0 = 1.1
#-----------------------------------------------------------------
#r = 0.58
#tau = np.linspace(0,math.pi,num = 500)
#t = np.linspace(0,150,300)
#k = np.sin(tau)+6.15

#plt.plot(tau,k)
#plt.xlabel('t')
#plt.ylabel('k')
#plt.title('Sinusoidal k(t) = sin(t)+6.15')
#plt.show()

#solvevaryk(g,t,tau,r,k)
#solvevaryk3(t,tau,r)
