import scipy.io 
import matplotlib as mpl
from scipy.optimize import fsolve
# Import math Library
import math 
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, solve_ivp
from sympy import Matrix
import random

def jordanFormsMat(A):
	k = A.shape[0]
	P = np.zeros((k,2,2), dtype = np.complex_)
	J = np.zeros((k,2,2), dtype = np.complex_)
	bigEigval = np.zeros ((k,1), dtype= np.complex_)
	for i in range(k):
		Ptemp, Jtemp=Matrix(A[i,:,:]).jordan_form()
		P[i,:,:] = np.array(Ptemp, dtype = np.complex_)
		#print(P[i,:,:])
		#Jtemp = np.real(Jtemp)
		J[i,:,:] = np.real(np.array(Jtemp, dtype = np.complex_))
		bigEigval[i] = max([J[i,0,0],J[i,1,1]])

		#print(J[i,:,:])
		#print(bigEigval[i])
		for j in range(A.shape[2]):
			colnorm = np.linalg.norm(P[i,:,j])
			P[i,:,j] =  P[i,:,j]/colnorm

	return P,J, bigEigval


def muG(P,J, G_adj):
	k = G_adj.shape[0]
	mu_G = float('-inf')
	sumln = 1
	sumnegl = 0
	sumposl = 0
	for i in range(k):
		for j in range(k):
			if G_adj[i,j]>0:
				prod = np.matmul(np.linalg.inv(P[j,:,:]),P[i,:,:])
				#print(P[j,:,:])
				#print(P[i,:,:])

				dist = np.linalg.norm(prod, ord=2)
				lmda = max([J[i,0,0],J[i,1,1]])
				print(lmda)
				if lmda <0:
					lmda = -lmda
					sumnegl = sumnegl + lmda
				else:
					sumposl = sumposl+lmda

				newmu_G = np.log(dist)/lmda
				sumln = sumln*dist

				if(newmu_G>mu_G):
					mu_G = newmu_G

	sumln = np.log(sumln)
	print(sumln)
	print(sumnegl)
	print(sumposl)
	return mu_G


def stabBoundRho2(P,J,G_adj):
	k = G_adj.shape[0]
	rho_2 = float('-inf')
	for i in range(k):
		for j in range(k):
			if G_adj[i,j]>0:
				prod = np.matmul(np.linalg.inv(P[j,:,:]),P[i,:,:])
				prod2 = np.matmul(np.linalg.inv(P[i,:,:]),P[j,:,:])

				dist = np.linalg.norm(prod, ord=2)
				dist2 = np.linalg.norm(prod2, ord = 2)

				lmdai = -max([J[i,0,0],J[i,1,1]])
				lmbdaj = -max([J[j,0,0],J[j,1,1]])

				newrho_2 = (np.log(dist)+np.log(dist2))/(lmdai+lmbdaj)
				print(newrho_2)
				if(newrho_2>rho_2):
					rho_2 = newrho_2
	return rho_2

def stabBoundRho(P,J,G_adj, G_cycles):
	n_cycles = G_cycles.shape[0]
	rho = float('-inf')
	for i in range(n_cycles):
		len_cycles = len(G_cycles[i])
		num = 0
		denom = 0
		for j in range(len_cycles-1):
			s = G_cycles[i][j+1]
			r = G_cycles[i][j]

			if G_adj[r,s]>0:
				#print(r,s)
				prod = np.matmul(np.linalg.inv(P[s,:,:]),P[r,:,:])

				dist = np.linalg.norm(prod, ord=2)

				lmdai = -max([J[r,0,0],J[r,1,1]])
				num = num + np.log(dist)
				denom = denom + lmdai

		new_rho = num/denom
		if(new_rho>rho):
				rho = new_rho
	return rho

def stabsldwell(P,J,G_adj,G_cycles):
	n_cycles = G_cycles.shape[0]
	nu = float('-inf')*np.ones((n_cycles,1))
	for i in range(n_cycles):
		len_cycles = len(G_cycles[i])
		num = 0
		denom = float('-inf')
		for j in range(len_cycles-1):
			s = G_cycles[i][j+1]
			r = G_cycles[i][j]

			prod = np.matmul(np.linalg.inv(P[s,:,:]),P[r,:,:])

			dist = np.linalg.norm(prod, ord=2)

			lmdai = np.real(max([J[r,0,0],J[r,1,1]]))
			num = num + np.log(dist)
			if(lmdai>denom):
				denom = lmdai

		nu[i] = num/-denom

	return nu

def sigmaSignal(N,dt,A,G_adj, bigEigval,tau, eta):
	sigma = np.zeros((N,1), dtype = int)
	dt_sigma = np.zeros((N,1))
	k = G_adj.shape[0]
	possiblenodes = range(k)
	currentnode = random.choice(possiblenodes)
	sigma[0] = int(currentnode)

	for i in range(N-1):
		if bigEigval[currentnode] >= 0:
			dt_sigma[i+1] = random.uniform(dt, eta)

		else:
			dt_sigma[i+1] = random.uniform(tau,2.5*tau)
		possiblenodes = np.where(G_adj[sigma[i],:][0]==1)
		possiblenodes = possiblenodes[0]
		currentnode = random.choice(possiblenodes)
		sigma[i+1] = int(currentnode)


	sigma = sigma[0: -2]
	dt_sigma = dt_sigma[1:-1]

	sigma = np.reshape(sigma,-1)
	dt_sigma = np.reshape(dt_sigma,-1)
	plt.plot(sigma,'o')
	plt.show()
	plt.plot(dt_sigma,'o')
	plt.show()
	t_sigma = np.cumsum(dt_sigma)
	T = t_sigma[-1]
	Nt_sigma = dt_sigma//dt
	Nt_sigma = Nt_sigma.astype(int)

	sigma_signal = np.repeat(sigma,Nt_sigma)
	t_signal = np.linspace(0,T,sigma_signal.size)
	plt.plot(t_signal, sigma_signal)
	plt.xlabel(r'time t')
	plt.ylabel(r'\sigma')
	plt.title(r'Signal \sigma')
	plt.show()
	return sigma, t_signal, dt_sigma


def f(x, *args):
	A = args[1]
	return np.matmul(A,x)
def trajsolve(sigma, t_signal, dt_sigma, dt, A, x0):
	y0 = x0
	last_t = 0
	#n_test = 300//dt
	#n_test = int(n_test)
	#t_signal = np.linspace(0, 300, n_test)
	fig, axs = plt.subplots(2)
	fig.suptitle('Trajectory of x(t)')
	for i in range(len(dt_sigma)-1):
		A_i = A[sigma[i]]
		tf = dt_sigma[i]
		N_i = tf//dt
		N_i = int(N_i)
		t = np.linspace(0,tf, N_i)

		sol = odeint(f,y0,t,args = (A_i,))
		axs[0].plot(t_signal[last_t:last_t+N_i],sol[:,0])
		axs[0].set(ylabel='x_1(t)')
		axs[1].plot(t_signal[last_t:last_t+N_i],sol[:,1])
		axs[1].set(xlabel='time t', ylabel='x_2(t)')

		#plt.plot(t_signal[last_t:last_t+N_i],sol[0,:])

		y0 = sol[-1]
		last_t = last_t+N_i-1
	plt.show()


A0 = np.array([[-1.5,0],[0,-1.5]])
A1 = np.array([[-1,0],[1,-1]])
A2 = np.array([[-11,3],[-18,4]])
A3 = np.array([[3,-45],[1,-11]])
A4 = np.array([[3,-46],[1,-11]])
A5 = np.array([[-2.1,1],[0,-2.1]])

A0 = np.array([[-1.5,0],[0,-1.5]])
A1 = np.array([[-1,0],[0,-1.1]])
A2 = np.array([[-11,3],[-18,4]])
A3 = np.array([[3,-45],[1,-11]])
A4 = np.array([[3,-46],[1,-11]])
A5 = np.array([[-2.1,0],[0,-2.2]])


A = np.array([A0,A1,A2,A3,A4,A5])

G_adj = np.array([[0, 1,0,0,0,0],[0,0,1,0,0,0],[1,0,0,1,1,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,1,0,0,0]])
G_cycles = np.array([[0,1,2,0],[2,3,2],[2,4,5,2]])

P,J, bigEigval = jordanFormsMat(A)

B0 = np.array([[-2,0],[0,-2]])
B1 = np.array([[-0.4,-0.03],[1.4,0.04]])
B2 = np.array([[0.9,0],[2,-0.6]])
B = np.array([B0,B1,B2])
G_adjB = np.array([[0, 1,0],[0,0,1],[1,0,0]])


PB,JB, bigEigvalB = jordanFormsMat(B)


#mu_G = muG(PB,JB,G_adjB)
#rho_2 = stabBoundRho2(P,J,G_adj)
rho = stabBoundRho(P,J,G_adj, G_cycles)
rho = np.real(rho)
#nu_C = stabsldwell(P,J,G_adj, G_cycles)
#print(rho_2)
#print(rho)
#print(nu_C)

N = 25
dt = 0.1
eta = 1
tau = rho
#sigma, t_signal , dt_sigma= sigmaSignal(N, dt,A, G_adj, bigEigval, tau, eta)

x0 = [0.1,0.2]
#trajsolve(sigma, t_signal, dt_sigma, dt, A, x0)
N = 25
tau = 2
eta = 0.8
sigma, t_signal , dt_sigma= sigmaSignal(N, dt,B, G_adjB, bigEigvalB, tau, eta)
trajsolve(sigma, t_signal, dt_sigma, dt, B, x0)
