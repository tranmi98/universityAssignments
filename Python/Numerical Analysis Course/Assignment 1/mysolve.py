# The function mysolve(A, b) is invoked by ccore.py
# to solve the linear system
# Implement your solver in this file and then run:
# python ccore.py -clscale float
# where the line argument 'clscale' allows the global coarsening of the mesh
# -clscale 1  fine mesh
# -clscale 10 coarse mesh

import numpy
import matplotlib
import time
# Return bool succes, sol x

def QR(A):
    size = numpy.shape(A);
    m = size[0];
    n = size[1];
    V = numpy.zeros((m, n), dtype = numpy.float64)
   
    for k in range(n):
        x = A[k:m, k]
        sign = numpy.sign(x[0])
        if sign == 0:
            sign = 1
        
        vk = x.copy()
        vk[0] = sign*numpy.linalg.norm(x) + x[0]
        vk = vk/numpy.linalg.norm(vk)
        A[k:m,k:n] = A[k:m, k:n] - 2*numpy.outer(numpy.matrix.transpose(numpy.array([vk])), numpy.array([numpy.dot(vk,A[k:m, k:n])]))
                
        V[k:m,k] = vk;
        

    return V, A
      

def QRsolve(R,b):
    size = numpy.shape(R);
    n = size[1]

   
    x = numpy.zeros((n,1))
    
    x[n-1] = b[n-1]/R[n-1,n-1]
   
    for i in range(n-2, -1, -1):
        term = numpy.dot(R[i, i+1 :n],x[i+1:n])
        x[i] = (b[i] - term)/R[i, i]
    return x  
def mysolve(A, b):
    
   
    A = numpy.array(A, dtype = numpy.float64)
    
    b = numpy.array(b, dtype = numpy.float64)
  
    V, R= QR(A)

    size = numpy.shape(R);
    m = size[0];
    n = size[1];
   
    for k in range(n):
        factor = numpy.dot(numpy.matrix.conjugate(V[k:m, k]),b[k:m])
        
        b[k:m] = b[k:m] - 2*numpy.array([V[k:m,k]])*factor
    
    x = QRsolve(R,b)

    return True, x
    
def genMat(n):
    M = numpy.random.rand(n,n)*100
    b = numpy.random.rand(n)*100
    if (n<200) and (numpy.linalg.det(M)==0):
        genMat(n)
    return M, b
def complexTemp():
    n = 8
    p = n/30
   
    axex = (2**numpy.arange(1,n,1))*10
    
    axey = numpy.zeros(n-1)
    compy = numpy.zeros(n-1)
    compx = (2**numpy.arange(1,n,1))*10

    q = 10
    for i in range(1,n,1):
        q = q*2

        M, b = genMat(q)
        before = time.time()
        bool, x = mysolve(M, b)
        after = time.time()
        axey[i-1] = after-before;
        if(i != 1):
            compy[i-1] = numpy.log2(axey[i-1]/axey[i-2])
        
        
        
        
    
    
    
    matplotlib.pyplot.plot(compx,compy)
    #matplotlib.pyplot.xscale('log', basex = 2)   
    #matplotlib.pyplot.yscale('log', basey = 2)  
    matplotlib.pyplot.xlabel("Taille de la matrice 2n")
    matplotlib.pyplot.ylabel("Complexité temporelle log2(O(2n)/O(n))")   
    matplotlib.pyplot.title("Complexité temporelle de resolution du systeme Ax = b par Householder")   
    matplotlib.pyplot.savefig("outComp2.png")

# decommenter pour executer les fonctions sur les complexites temporelles
#complexTemp()
