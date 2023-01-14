
"""
Created on Sat Dec 22 14:58:26 2018

@author: minh-phuong
"""

import numpy as np
import ILU0
def csrProd(sA, iA, jA, p):
    n = len(p)
    result = np.zeros(n)
    for i in range(n):
        for j in range(iA[i], iA[i+1]):
            result[i] = result[i] + sA[j]*p[jA[j]]
    return result

def csrCG (sA,iA,jA,b,rtol,prec):
    
    n = len(b)
    x = np.zeros(n)
    r = np.copy(b)
    p = np.copy(r)
    #preallouer taille n+1 car convergence <= n
    res = np.zeros(n+1)
    
    k = 0
    res[k] = np.linalg.norm(r)
    if prec == False:
        while res[k]>rtol:
            prod = csrProd(sA, iA, jA, p)
            alpha = np.dot(r, r)/np.dot(p, prod)
            x = x + alpha*p
            oldr = np.copy(r)
            r = r - alpha*prod
            beta = np.dot(r, r)/np.dot(oldr, oldr)
            p = r + beta*p
            k = k+1
            res[k] = np.linalg.norm(r)
    else:
        sLU, iLU, jLU = ILU0.LUcsr(np.copy(sA),np.copy(iA),np.copy(jA))
        y = ILU0.fwdSubcsr(sLU, iLU, jLU, r)
        rtilde = np.reshape(ILU0.bwdSubcsr(sLU, iLU, jLU, y), len(y))
        p = np.copy(rtilde)
        while res[k]>rtol:
            prod = csrProd(sA, iA, jA, p)
            alpha = np.dot(r, rtilde)/np.dot(p, prod)
            x = x + alpha*p
            oldr = np.copy(r)
            oldrtilde = np.copy(rtilde)
            r = r - alpha*prod
            y = ILU0.fwdSubcsr(sLU, iLU, jLU, r)
            rtilde = np.reshape(ILU0.bwdSubcsr(sLU, iLU, jLU, y), len(y))
            beta = np.dot(r, rtilde)/np.dot(oldr, oldrtilde)
            p = rtilde + beta*p
            k = k+1
            if (k == len(res)):
                new = np.zeros(k*2)
                new[:k] = res
                res = new
            res[k] = np.linalg.norm(r)
            
    res = res[:k+1]
    return (x, res)
