import numpy as np
import sys
from math import pi

# Number of grid points and mesh size
n=100
h=1/n
f=1/(h*h)

# Nonlinear function in the BVP
def fu(z):
    return 80*np.cos(z)

# Derivative of the nonlinear function, used to compute the Jacobian
def dfu(z):
    return -80*np.sin(z)

# Initial guess for the solution
U=np.zeros((n-1))

# Perform Newton iterations
F=np.empty((n-1))
tol=1e-10
k=0
while k==0 or np.linalg.norm(dU)/np.linalg.norm(U)>tol:

    # Compute the Jacobian
    A=np.diag(np.ones(n-1)*(-f*2-dfu(U)))+np.diag(np.ones(n-2)*f,1)+np.diag(np.ones(n-2)*f,-1)

    # Compute F(U)
    for i in range(0,n-1):
        F[i]=-2*f*U[i]-fu(U[i])
        if i>0: F[i]+=f*U[i-1]
        if i<n-2: F[i]+=f*U[i+1]
        else: F[i]+=f*10

    # Compute the Newton update
    dU=np.linalg.solve(A,F);
    U-=dU
    k+=1

# Output the solution, appending terms at the x=0 and x=1 that were not solved
# for in the linear system
print(0,0)
for i in range(n-1):
    print((i+1)*h,U[i])
print(1,10)
