
"""
Created on Sat Dec 22 15:00:24 2018

@author: minh-phuong
"""
import numpy as np
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

def LUcsr(sA, iA, jA):
    sLU = sA
    iLU = np.array(iA).astype(int)
    jLU = np.array(jA).astype(int)
    iLU = np.reshape(iLU, len(iLU))
    jLU = np.reshape(jLU, len(jLU))
    n=len(iA)-1
    for i in range(1, n):
        rowi = np.array(range (iA[i],iA[i+1]))
        rowik = rowi[jA[rowi]<i]
        for indk in rowik:
            columnk = jA[indk]
            rowk = np.array(range(iA[columnk], iA[columnk+1]))
            skk = rowk[jA[rowk]== columnk]
            sLU[indk] = sLU[indk]/sLU[skk]
            for indj in range(indk+1, iA[i+1]):
                columnj = jA[indj]
                skj = rowk[jA[rowk]== columnj]
                if len(skj)!=0:
                    sLU[indj] = sLU[indj] - sLU[indk]*sLU[skj]
    
    return sLU, iLU, jLU