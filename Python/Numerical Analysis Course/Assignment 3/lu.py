# The function mysolve(A, b) is invoked by ccore.py
# to solve the linear system
# Implement your solver in this file and then run:
# python ccore.py -clscale float
# where the line argument 'clscale' allows the global coarsening of the mesh
# -clscale 1  fine mesh
# -clscale 10 coarse mesh


SolverType = 'scipy'

import numpy as np
import matplotlib as pl
import time 

lowband = 0
highband = 0
#########################################################################################
###### Outils
#########################################################################################    
def density(A):
    n = A.shape[0]
    return sum(sum(A!=0))/(n*n)

def densitycsr(sA,iA):
    n = len(iA)-1
    
    return len(sA)/(n*n)


#########################################################################################
###### Matrice Pleine
#########################################################################################
def LUfactorize(A):
   
   m = len(A)   
   for k in range(m-1): 
       alpha = A[k,k]
       for j in range(k+1,m):
           A[j,k] = A[j,k]/alpha
           A[j,k+1:m] = np.subtract(A[j,k+1:m],np.dot(A[j,k],A[k,k+1:m]))



def fwdSub(A,b):
    
    m = len(b)
    y = np.array(np.zeros((m,1)), dtype = np.float64)
    y[0] = b[0]
    for i in range(1,m):
        term = np.dot(A[i, 0:i], y[0:i])
        y[i] = b[i] - term
        
    return y

def bwdSub(A,b):
    n = len(b);
    x = np.zeros((n,1))
    x[n-1] = b[n-1]/A[n-1,n-1]  
    for i in range(n-2, -1, -1):
        term = np.dot(A[i, i+1 :n],x[i+1:n])
        x[i] = (b[i] - term)/A[i, i]
        
    return x  

def LUsolve(A, b):
    A = np.array(A, dtype = np.float64)        
    b = np.array(b, dtype = np.float64)
    before = time.time()

    LUfactorize(A) 
    denLU = density(A)
    y = fwdSub(A,b)
    x = bwdSub(A,y)
    after = time.time()
    return np.transpose(x)
##########################################################################################
###### CSR
##########################################################################################
    
def bandCSR(iA, jA):
    n = len(iA)-1
    high = jA[iA[1::]-1] - range(n)
    low = jA[iA[:-1:]] - range(n)
    highband = max(high)
    lowband = min(low)
    return lowband, highband

def CSRformat(A):

    [m,n] = np.shape(A)

    indnonzero = np.nonzero(A)
    sA = A[indnonzero]
    jA = indnonzero[1]
    rows = indnonzero[0]
    iA = [0]*(m+1)
    iA = np.array(iA)
    for i in rows:
        iA[i+1::] = iA[i+1::]+1

    iA = np.array(iA).astype(int)
    jA = np.array(jA).astype(int)
    return sA, iA, jA

def CSRtoBand(sA, iA, jA):
    global lowband, highband
    lowband, highband = bandCSR(iA, jA)
    m = len(iA)-1
    bandwidth = highband - lowband + 1 
    sAp = np.array(np.zeros((m*bandwidth,1)), dtype = np.float64)
    jAp = np.zeros((m*bandwidth,1))
    iAp = np.zeros((m+1,1))
    count = 0
    for i in range(m):

        left = max(i+lowband, 0)
        right = min(i+highband+1, m)
        bandi = right - left
        for j in range(iA[i],iA[i+1]):
            sAp[count+jA[j]-left] = sA[j]
        
        col = np.array(range(left,bandi+left))
        jAp[count:count+bandi] = np.reshape(col,(len(col), 1))
        count = count + bandi
        iAp[i+1] = count

    
    sAp = np.delete(sAp, np.s_[count:])
    jAp = np.delete(jAp, np.s_[count:])
    return sAp, iAp, jAp

def nonNul(sA,iA,jA):
    sLU = sA
    iLU = iA.astype(int)
    jLU = jA.astype(int)
    
    indnonzero = np.nonzero(sLU)
    indnonzero = indnonzero[0]
    sLU = sLU[indnonzero]
    jLU = jLU[indnonzero]
        
    one = np.ones(len(iA)).astype(int)
    term = np.zeros(len(iA)).astype(int)
    notin = np.array(range(len(sA)))
    notin = np.delete(notin, indnonzero)
    for i in notin:
        ind = jA[i].astype(int)
        term[ind+1::] = term[ind+1::]+ one[ind+1::]
    
    term = np.reshape(term, np.shape(iLU))

    iLU = iLU - term
    return sLU, iLU, jLU

def LUcsr(sAp, iAp, jAp):
    sA,iA,jA = CSRtoBand(sAp,iAp,jAp)
    global lowband, highband
    sLU = sA
    iLU = np.array(iA).astype(int)
    jLU = np.array(jA).astype(int)
    iLU = np.reshape(iLU, len(iLU))
    jLU = np.reshape(jLU, len(jLU))

    m = len(iA)-1

    for k in range(m-1):
        indiag = min(k, abs(lowband))
        lowk = iLU[k]+indiag
        alpha = sLU[lowk]
        
        down = min(k-lowband +1 , m)
        for j in range(k+1,down):

            indiagj = min(j, abs(lowband))
            dif = j - k

            low = iLU[j]+indiagj - dif
            sLU[low] = sLU[low]/alpha

            high = min(low + 1 + highband, iLU[j+1])
            highk = min(lowk+ 1+ highband, iLU[k+1])


            term = np.dot(sLU[low],sLU[lowk + 1: highk])
            sLU[low + 1 :high] = np.subtract(sLU[low + 1 :high],term)
    
    
    sLU, iLU, jLU = nonNul(sLU, iLU, jLU)     
   # sLU1,iLU1,jLU1 = nonNul1(sLU, iLU, jLU)
   # print(np.all(sLU1==sLU2),np.all(jLU1==jLU2), np.all(iLU1==iLU2))
    return sLU, iLU, jLU 
        
def fwdSubcsr(sA, iA, jA, b):

    m = len(b)
    y = np.array(np.zeros((m,1)), dtype = np.float64)
    sA = np.reshape(sA, (len(sA),1))
    y[0] = b[0]
    for i in range(1,m):
        start = iA[i]
        fin = iA[i+1]
        line = jA[start:fin]
        end = np.where(line==i)
        end = start+end[0][0]
        length = end-start
        prod = 0
        if(length != 0):
            xind = jA[start:end]
            prod = np.dot(np.transpose(sA[start:end]), y[xind])
        
        y[i] = b[i] - prod
        
    return y                

def bwdSubcsr(sA, iA, jA, b):
    n = len(b);
    x = np.array(np.zeros((n,1)), dtype = np.float64)
    x[n-1] = b[n-1]/sA[-1]  
    for i in range(n-2, -1, -1): 
        start = iA[i]
        fin = iA[i+1]
        line = jA[start:fin]
        begin = np.where(line==i)
        begin = start+begin[0][0]
        length = fin - begin
        alpha = sA[begin]
        prod = 0
        if(length !=0):
            xind= jA[begin+1:fin]
            prod = np.dot(np.transpose(sA[begin+1:fin]), x[xind])
        
        x[i] = (b[i] - prod)/alpha
        
    return x  

def LUcsrsolve(A, b):
    
    A = np.array(A, dtype = np.float64)

    b = np.array(b, dtype = np.float64)
    sA, iA, jA = CSRformat(A)
    before = time.time()
    sLU, iLU, jLU = LUcsr(sA,iA,jA)
    denLU = densitycsr(sLU, iLU)
    y = fwdSubcsr(sLU, iLU, jLU, b)
    x = bwdSubcsr(sLU, iLU, jLU, y)
    after = time.time()
    return np.transpose(x)
##########################################################################################
##### RCMK
###########################################################################################
def RCMK(iA,jA):
    n = len(iA)-1
    iA = np.reshape(iA.astype(int), len(iA))
    jA = np.reshape(jA.astype(int), len(jA))
    r = np.array([None]*n)
    waiting = [None]*n
    waiting = np.array(waiting)
    totake = [None]*(n+1)
    totake[:n:] = range(n)
    totake = np.array(totake)
    degrees = np.subtract(iA[1::],iA[:-1:])
    degrees = np.reshape(degrees.astype(int), len(degrees))
     
    countwait = 0
    countr = 0

   
    while not np.all(waiting==None) or not np.all(totake==None):

        if not np.all(waiting==None):
            node = waiting[0].astype(int)
            waiting[0] = None
            waiting = np.roll(waiting, -1)
            countwait = countwait - 1
            r[countr] = node
            countr = countr + 1
            numneigh = degrees[node]
            neighbors= jA[range(iA[node], iA[node+1])]
            neighi = [None]*(numneigh+1)
            neighi = np.array(neighi)
            neighi[:numneigh:] = iA[neighbors+1] - iA[neighbors]
            
            sortneigh = np.sort(neighi[:numneigh:])
            if(numneigh >1):
                for i in range(numneigh):
                   nextint = np.where(neighi==sortneigh[i])
                   nextnode = jA[iA[node]+nextint[0][0]]
                   neighi[nextint[0][0]] = None

                   if(nextnode not in r and nextnode not in waiting):
                       waiting[countwait] = nextnode
                       totake[nextnode] = None
                       countwait = countwait + 1                

        elif not np.all(totake==None):
            notnone = np.array(np.where(totake != None))
            notnone = np.reshape(notnone, notnone.size)
            newdegrees = np.subtract(iA[notnone+1],iA[notnone])
            mindeg = notnone[np.argmin(newdegrees)]
            waiting[countwait] = mindeg
            countwait = countwait + 1
            totake[mindeg] = None


    r = r[::-1]    
    return r

###AVEC L'AIDE DE ANTOINE TUERLINCKX
def RCMKchange(sA, iA, jA, r):
    
    n=len(iA)-1
    iAp = np.zeros(n+1).astype(int)
    sAp=np.zeros(len(sA))
    jAp=np.zeros(len(sA)).astype(int)
    
    rind = np.array(range(len(r)))
    rinv = np.zeros(len(r)).astype(int)
    rinv[r] = rind
    
    difi = np.array(iA)
    difi = difi[1::] - difi[:-1:]
    difi= difi[r]
    
    for i in range(n):

        start = iA[r[i]]
        end = iA[r[i]+1]
        
        iAp[i+1] = iAp[i]+difi[i]
        
        columns=rinv[jA[start:end]]
        sAline=sA[start:end]
        zipped = zip(columns, sAline)
        sort = sorted(zipped)
        columns,sAline = zip(*sort)
        sAline = np.array(sAline)
        columns = np.array(columns)
        
        sAp[iAp[i]:iAp[i+1]]=sAline
        jAp[iAp[i]:iAp[i+1]]=columns
    
    
    return sAp,iAp,jAp

def LUcsrRCMKsolve(A,b):
    
    A = np.array(A, dtype = np.float64)

    b = np.array(b, dtype = np.float64)
    sAp, iAp, jAp = CSRformat(A) 
    iAp = np.array(iAp)
    jAp = np.array(jAp)
    bef = time.time()
    r = RCMK(iAp,jAp)    
    r = np.array(r).astype(int)
    af = time.time()
    sA2, iA2, jA2 = RCMKchange(sAp,iAp,jAp, r)

    before = time.time()

    sLUd, iLUd, jLUd = LUcsr(sA2,iA2,jA2)
    denLU = densitycsr(sLUd, iLUd)
    bd = b[r]
    yd = fwdSubcsr(sLUd, iLUd, jLUd, bd)
    xd = bwdSubcsr(sLUd, iLUd, jLUd, yd)
    xt = np.array([0.0]*len(xd))
    r = np.reshape(r, (len(xd),1))
    xt[r] = xd
    after = time.time()    
    
    return np.transpose(xt)



