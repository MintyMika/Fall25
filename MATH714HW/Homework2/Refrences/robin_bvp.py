#!/usr/bin/python3
import numpy as np
from math import *

# Reference solution
def uex(z):
    return -0.5*exp(pi)*(sin(z)+cos(z))-0.5*exp(z)

# Solves the BVP with the Robin boundary conditions
def robin_bvp(n,show):

    # Set up grid and source term
    h=pi/n
    x=np.linspace(0,pi,n+1)
    F=-h*h*np.exp(x)

    # Generate the centered difference differentiation matrix for u''
    A=np.diag(np.ones(n+1)*(h*h-2))+np.diag(np.ones(n),1)+np.diag(np.ones(n),-1)

    # Modify the first and last rows of the matrix to account for the Robin
    # boundary conditions
    A[0,0]=h*h-2-2*h
    A[0,1]=2
    A[n,n]=h*h-2-2*h
    A[n,n-1]=2

    # Solve the linear system
    U=np.linalg.solve(A,F)

    # Optionally print out the solution
    if(show):
        for j in range(n+1):
            uex_=uex(j*h)
            print(x[j],U[j],uex_,U[j]-uex_)

    # Compute the 2-norm error
    err=0.5*((U[0]-uex(0))**2+(U[n]-uex(pi))**2)
    for j in range(1,n):
        err+=(U[j]-uex(h*j))**2
    return sqrt(h*err)

# Optional line to print out the numerical solution in comparison to the
# reference solution 
#robin_bvp(80,True)

# Calculate the 2-norm errors for a range of grid sizes
e20=robin_bvp(20,False)
print(20,e20,1)
for n in (40,80,160):
    err=robin_bvp(n,False)
    print(n,err,e20/err)
