# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 09:04:10 2018

@author: sebco and minh-phuong tran
"""

import scipy.linalg
import numpy  as np
import matplotlib

import lu
import time
eps = 1.e-6

def complexTemp():
    n = 13
    p = 2
    q = 1

    axex = p*(np.arange(1,n,1))+q
    
    yLU = np.zeros(n-1)
    yCSR = np.zeros(n-1)
    yRCMK = np.zeros(n-1)
    transCSR = np.zeros(n-1)
    
    compLU = np.zeros(n-1)
    compCSR = np.zeros(n-1)
    compRCMK = np.zeros(n-1)
    compformCSR = np.zeros(n-1)

    compx = p*(numpy.arange(1,n,1))+q

    for i in range(1,n,1):
        q = q+p
        m = q*q
        A = creatDF(q)
        b = 5*np.random.rand(m)-10

        before = time.time()
        xtest = lu.LUsolve(np.copy(A),b)
        after = time.time()
        yLU[i-1] = after-before;
        
        before = time.time()
        sA,iA,jA = lu.CSRformat(np.copy(A))
        after = time.time()
        transCSR[i-1] = after-before;
        
        before = time.time()
        xtest = lu.LUcsrsolve(np.copy(A),b)
        after = time.time()
        yCSR[i-1] = after-before;
        
        before = time.time()
        xtest = lu.LUcsrRCMKsolve(np.copy(A),b)
        after = time.time()
        yRCMK[i-1] = after-before;
        
        q2 = q*2
        m2 = q2*q2
        
        A2 = creatDF(q2)
        b2 = 5*np.random.rand(m2)-10

        before = time.time()
        xtest = lu.LUsolve(np.copy(A2),b2)
        after = time.time()
        time2LU = after-before;
        
        before = time.time()
        xtest = lu.LUcsrsolve(np.copy(A2),b2)
        after = time.time()
        time2CSR = after-before;
        
        before = time.time()
        xtest = lu.LUcsrRCMKsolve(np.copy(A2),b2)
        after = time.time()
        time2RCMK = after-before;
        
        before = time.time()
        sA,iA,jA = lu.CSRformat(np.copy(A2))
        after = time.time()
        time2 = after-before;
        
        if(i != 1):
            compLU[i-1] = np.log(time2LU/yLU[i-1])/np.log(4)
            compCSR[i-1] = np.log(time2CSR/yCSR[i-1])/np.log(4)
            compRCMK[i-1] = np.log(time2RCMK/yRCMK[i-1])/np.log(4)
            compformCSR[i-1] = np.log(time2/transCSR[i-1])/np.log(4)

        
        
        
    
    
    matplotlib.pyplot.figure(1)
    matplotlib.pyplot.plot(compx,compLU, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Complexité temporelle log4(O(4n)/O(n))")   
    matplotlib.pyplot.title("Complexité temporelle de resolution du systeme Ax = b par LU plein")   
    matplotlib.pyplot.savefig("compLU4.png")

    matplotlib.pyplot.figure(2)
    matplotlib.pyplot.plot(compx,compCSR, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Complexité temporelle log4(O(4n)/O(n))")   
    matplotlib.pyplot.title("Complexité temporelle de resolution du systeme Ax = b par LU creux")   
    matplotlib.pyplot.savefig("compCSR4.png")
    
    matplotlib.pyplot.figure(3)
    matplotlib.pyplot.plot(compx,compLU, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Complexité temporelle log4(O(4n)/O(n))")   
    matplotlib.pyplot.title("Complexité temporelle de resolution du systeme Ax = b par LU creux avec RCMK")   
    matplotlib.pyplot.savefig("compRCMK4.png")

    matplotlib.pyplot.figure(4)
    matplotlib.pyplot.plot(axex,yLU, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Temps de résolution")   
    matplotlib.pyplot.title("Temps de resolution du systeme Ax = b par LU plein")   
    matplotlib.pyplot.savefig("tLU4.png")

    matplotlib.pyplot.figure(5)
    matplotlib.pyplot.plot(axex,yCSR, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Temps de résolution")   
    matplotlib.pyplot.title("Temps de resolution du systeme Ax = b par LU creux")   
    matplotlib.pyplot.savefig("tCSR4.png")
    
    matplotlib.pyplot.figure(6)
    matplotlib.pyplot.plot(axex,yRCMK, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Temps de résolution")   
    matplotlib.pyplot.title("Temps de resolution du systeme Ax = b par LU creux avec RCMK")   
    matplotlib.pyplot.savefig("tRCMK4.png")
    
    matplotlib.pyplot.figure(7)
    matplotlib.pyplot.plot(axex,transCSR, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Temps d'exécution")   
    matplotlib.pyplot.title("Calcul du format CSR")   
    matplotlib.pyplot.savefig("csrform.png")
    
    matplotlib.pyplot.figure(8)
    matplotlib.pyplot.plot(axex,compformCSR, '-bo')
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n^2")
    matplotlib.pyplot.ylabel("Complexité temporelle log4(O(4n)/O(n))")   
    matplotlib.pyplot.title("Complexité temporelle du calcul du format CSR")  
    matplotlib.pyplot.savefig("csrformcomp.png")


def complexTempLUsolve():
    n = 12
    f = 5
    axex = f*(np.arange(1,n))
    
    axey = np.zeros(n-1)
    compy = np.zeros(n-1)
    compx = f*(np.arange(1,n))
    j = 1
    for i in range(1,n):
        j = 5*i
        m = j*j
        A = creatDF(j)
        b = 5*np.random.rand(m)-10
        before = time.time()
        xtest = lu.LUcsrRCMKsolve(np.copy(A),b)
        after = time.time()
  
        axey[i-1] = after-before;
        if(i != 1):
            compy[i-1] =np.log(axey[i-1]/axey[i-2]) / np.log(4)
        
        
        
        
    
    
    
    matplotlib.pyplot.plot(axex,axey)
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice n")
    matplotlib.pyplot.ylabel("Temps de résolution")   
    matplotlib.pyplot.title("Temps de résolution avec décomposition LU creuse avec RCMK")   
    matplotlib.pyplot.savefig("LUrcmktemp.png")
    
"""TEST LUfactorize"""
def testLUfactorize():
    print('Test LUfactorize')
    print('---------------------------------------')

    #Matrix
    n=100
    A = creatDF(10)
    
    #Scipy LU
    [P,L,U] = scipy.linalg.lu(np.copy(A))
    
    #your LU
    LU = np.copy(A)
    lu.LUfactorize(LU)
    
    Ltest = np.eye(n) + np.tril(LU,-1)
    Utest = np.triu(LU)
    
    Atest = Ltest @ Utest
    if(np.all(abs(Atest - A)<= eps)):
        print('LU decomposition is exact')
    else :
        print('Mistake found in LU decomposition')
    if(np.all(abs(Ltest - L)<= eps)):
        print('L is exact')
    else :
        #print(abs(Ltest - L))
        print('L not exact')
        
    if(np.all(abs(Utest - U)<= eps)):
        print('U is exact')
    else :
        #print(abs(Utest - U))
        print('U not exact')
    print('---------------------------------------\n')
    
"""Test LUsolve"""
def testLUsolve():
    print('Test LUsolve')
    print('---------------------------------------')
    
    n = 100
    A = creatDF(10)
    b = 5*np.random.rand(n)-10
    x = scipy.linalg.solve(np.copy(A),b)
    before = time.time()
    xtest = lu.LUsolve(np.copy(A),b)
    after = time.time()
    print(after - before)
    if(np.all( abs(xtest - x) <= eps)):
        print('Solver works correctly : solution x is exact')
    else :
        print('Solver does not work correctly : solution x is not exact')
    print('---------------------------------------\n')

"""TEST CSR Format"""
def testCSRformat():
    print('Test CSRformat')
    print('---------------------------------------')
    #1st Example
    A1 = np.zeros((4,4))
    [sA1,iA1,jA1] = lu.CSRformat(A1)
    iA1_sol = np.zeros(5,dtype = int)
    jA1_sol = np.zeros(0,dtype = int)
    sA1_sol = np.zeros(0,dtype = complex)
    if (np.all(iA1 == iA1_sol) and np.all(jA1 == jA1_sol) and np.all(sA1 == sA1_sol)):
        print('CSR Format is exact for example 1')
    else :
        print('CSR Format is not correct for example 1')
        print(iA1,jA1,sA1)
        print(iA1_sol,jA1_sol,sA1_sol)
    
    #2nd Example
    A2 = np.array([[1,2,0],[4,5,6],[0,0,1]])
    [sA2,iA2,jA2] = lu.CSRformat(A2)
    iA2_sol = np.array([0,2,5,6],dtype = int)
    jA2_sol = np.array([0,1,0,1,2,2],dtype = int)
    sA2_sol = np.array([1,2,4,5,6,1],dtype = complex)
    
    if (np.all(iA2 == iA2_sol) and np.all(jA2 == jA2_sol) and np.all(sA2 == sA2_sol)):
        print('CSR Format is exact for example 2')
    else :
        print('CSR Format is not correct for example 2')
        print(iA2,jA2,sA2)
        print(iA2_sol,jA2_sol,sA2_sol)
    print('---------------------------------------\n')

    
"""TEST LUCSR"""
def testLUcsr():
    print('Tesr LUcsr')
    print('---------------------------------------')
    
    #EX1 : small matrice - not sparse - not symmetric but positive definite
    A = np.array([[3.,-1.,1.],[-2.,2.,-1.],[0.,-1.,4.]])
    #print(A)
    #print(np.linalg.eigvals(A))

    #Your solution
    [sA,iA,jA] = lu.CSRformat(np.copy(A))
    [sLU,iLU,jLU] = lu.LUcsr(sA,iA,jA)
    
    #True solution
    [P,L,U] = scipy.linalg.lu(A)
    LU_sol = np.tril(L,-1) + U
    [sLU_sol,iLU_sol,jLU_sol] = lu.CSRformat(LU_sol)
    #print(LU_sol)

    if (np.all(iLU == iLU_sol) and np.all(jLU == jLU_sol) and np.all(abs(sLU - sLU_sol)<=eps)):
        print('LUcsr is exact for example 1 : small matrice (not sparse - not symmetric) ')
    else:
        print(sLU_sol)
        print(sLU)

        print('LUcsr is not correct for example 1 : small matrice (not sparse  - not symmetric)')
       
    # EX2 : bigger matrix but sparse, symetric and positive definite
    A = creatDF(10) #create matrix 100*100

    #Your solution
    [sA,iA,jA] = lu.CSRformat(np.copy(A))
    [sLU,iLU,jLU] = lu.LUcsr(sA,iA,jA)
    #True solution
    [P,L,U] = scipy.linalg.lu(A)
    LU_sol = np.tril(L,-1) + U
    [sLU_sol,iLU_sol,jLU_sol] = lu.CSRformat(LU_sol)
    
    if (np.all(iLU == iLU_sol) and np.all(jLU == jLU_sol) and np.all(abs(sLU - sLU_sol)<=eps)):
        print('LUcsr is exact for example 2 : bigger matrix but sparse and symmetric')
    else :
        print('LUcsr is not correct for example 2 : bigger matrix but sparse an symmetric')
        print(sLU_sol)
        print(sLU)
    print('---------------------------------------\n')

"""Test LUcsrsolve"""
def testLUcsrsolve():
    print('Test LUcsrsolve')
    print('---------------------------------------')
    n = 1600
    A = creatDF(40)
    #print(A)
    b = 5*np.random.rand(n)-10
    x = scipy.linalg.solve(np.copy(A),np.copy(b))
    
    before = time.time()
    xtest = lu.LUcsrsolve(np.copy(A),np.copy(b))
    after = time.time()
    print(after - before)
    

    #y = np.dot(A,x)
    #ytest = np.dot(A,np.reshape(xtest, len(xtest[0])))


    if(np.all( abs(xtest - x) <= 1.e-6)):
        print('Solver works correctly : solution x is exact')
    else :
        print('Solver does not work correctly : solution x is not exact')
        print(x)
        print(xtest)
        #print(y)
       # print(ytest)
    print('---------------------------------------\n')
    
    """Test LUcsrRCMKsolve"""
def testLUcsrRCMKsolve():
    print('Test LUcsrRCMKsolve')
    print('---------------------------------------')
    n = 16
    A = creatDF(4)
    #print(A)
    b = 5*np.random.rand(n)-10
    x = scipy.linalg.solve(np.copy(A),np.copy(b))
    
    before = time.time()
    xtest = lu.LUcsrRCMKsolve(np.copy(A),np.copy(b))
    after = time.time()
    print(after - before)
    

    #y = np.dot(A,x)
    #ytest = np.dot(A,np.reshape(xtest, len(xtest[0])))


    if(np.all( abs(xtest - x) <= 1.e-6)):
        print('Solver works correctly : solution x is exact')
    else :
        print('Solver does not work correctly : solution x is not exact')
        print(x)
        print(xtest)
        #print(y)
       # print(ytest)
    print('---------------------------------------\n')
    
"""Test 3LUsolve"""
def test3LUsolve():
    print('Test 3LUsolve')
    print('---------------------------------------')
    n = 1600
    A = creatDF(40)
    b = 5*np.random.rand(n)-10
    x = scipy.linalg.solve(np.copy(A),b)
    print("full")
    before = time.time()
    xtest = lu.LUsolve(np.copy(A),b)
    after = time.time()
    print("full time",after-before)
    print("csr")
    before = time.time()
    xtest1 = lu.LUcsrsolve(np.copy(A),np.copy(b))
    after = time.time()
    print("csr time",after - before)
    print("rcmk")
    before = time.time()
    xtest2 = lu.LUcsrRCMKsolve(np.copy(A),np.copy(b))
    after = time.time()
    print("rcmk time", after - before)
    if(np.all( abs(xtest - x) <= eps) and np.all( abs(xtest1 - x) <= eps) and np.all( abs(xtest2 - x) <= eps)):
        print('Solver works correctly : solution x is exact')
    else :
        print('Solver does not work correctly : solution x is not exact')
    print('---------------------------------------\n')    
"""
creatSymDF crée une matrice A de différences finies de taille n²xn² pour une problème de Laplace
@param : n taille de la discrétisation du problème de Laplace
@return : la partie symétrique de A, définie positive dimension n² x n²
"""
def creatDF(n):
    n2 = n * n
    A = np.zeros((n2,n2))
    for i in range(1,n+1):
        for j in range(1,n+1):
            index = i + (j-1)*n - 1

            if i==1 or i==n or j==1 or j==n :
                A[index,index] = 1.
            else :
                A[index,index] = 4.
                A[index,index+1] = -1.
                A[index,index-1] = -1.
                A[index,index+n] = -1.
                A[index,index-n] = -1.
    return A @ A.T


"""
density calcule la densité d'une matrice
@param : A est une matrice passée sous forme de numpy array
@return : densité de A (nnz/(n*n))
"""

def density(A):
    n = A.shape[0]
    return sum(sum(A!=0))/(n*n)


#TESTS
testLUfactorize()
testLUsolve()
testCSRformat()
testLUcsr()
testLUcsrsolve()
testLUcsrRCMKsolve()
test3LUsolve()
complexTemp()
