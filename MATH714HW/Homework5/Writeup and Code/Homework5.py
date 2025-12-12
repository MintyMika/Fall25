import numpy as np
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import quad



def Problem1():
    # set up the different stencils for easy access.
    stencil_coeffs=np.array([[3./8,-10./8,15./8,0,0],
                            [0,-1./8,6./8,3./8,0],
                            [0,0,3./8,6./8,-1./8]])

    beta_coeffs=np.array([[4./3,-19./3,25./3,11./3,-31./3,10./3],
                        [4./3,-13./3,13./3,5./3,-13./3,4./3],
                        [10./3,-31./3,25./3,11./3,-19./3,4./3]])
        
    # short routine for calculating the smoothness indicator bk.
    def smoothness(k,ul,uc,ur):
        c=beta_coeffs[k]
        us=np.array([ul**2,ul*uc,uc**2,ul*ur,uc*ur,ur**2])
        return np.dot(c,us)

    def eno(um2,um1,u,up1,up2):
        us = np.array([um2,um1,u,up1,up2])
        '''us - the set of 5 values at [j-2,j-1,j,j+1,j+2] of the density
        from which to construct the interpolated value at j+1/2.
        Complete the function eno, which should compute the smoothness
        indicator beta of each of the three possible stencils from these
        five points, choose the scheme with the lowest beta, and output
        the interpolated value u_j+1/2.
        '''
        # compute the smoothness indicator of each stencil.
        beta=np.zeros(3)
        for k in range(3):
            beta[k]=smoothness(k,us[k],us[k+1],us[k+2])
        opt=np.argmin(beta) # choose scheme with lowest beta
        c=stencil_coeffs[opt]
        return np.dot(c,us) # multiply density by the stencil coefficients
    
    # Grid size
    m=201
    # Pad by 6 ghost nodes: three on left, three on right. 
    # For c>0 only one set is used, but we can maintain a
    # more general implementation here. Also, we will take
    # advantage of Python's negative indexing and put all 6
    # ghost nodes at the end of the interval, so the indexing
    # across the array will go
    # 0,   1,   ...,   j,   ...,   m-1,   m,   m+1,   m+2,   -3,   -2,   -1.
    #|_________domain [-1,1]___________||___right ghost___||___left ghost __|
    mp6=m+6

    # PDE-related constants.
    c=1.
    dx=2.0/(m-1)
    dt=0.001
    nu=c*dt/dx
    T=1.5 # final simulation time
    snaps=3 # number of snapshots to output (excluding t=0)
    iters=int(T/dt)//snaps # iterations to perform between snapshots

    u=np.empty(mp6) # memory for the current step
    u1=np.empty(mp6) # memory for the next step
    us=np.empty((m,snaps+1)) # memory for all snapshots

    # Initial condition
    k=100
    for j in range(-3,m+3,1): # include ghost nodes
        x=-1+dx*j
        u[j]= 1./(1+np.exp(-k*x))

    u1=np.copy(u)
    us[:,0]=u[:m]


    
    # Integrate the PDE
    for i in range(1,snaps+1):
        for k in range(iters):
            for j in range(m):
                ur=eno(u[j-2],u[j-1],u[j],u[j+1],u[j+2]) # u_j+1/2
                ul=eno(u[j-3],u[j-2],u[j-1],u[j],u[j+1]) # u_j-1/2
                u1[j]=u[j]-nu*(ur-ul)
            u=np.copy(u1)
        us[:,i]=u[:m]



    # Run two cases: red light and green light
    def run_case(q_l, q_r, case_name):
        # Initialize solution
        for j in range(m):
            x = -1 + dx * j
            u[j] = q_r if x >= 0 else q_l
        
        # Set ghost nodes
        u[-3:] = q_l  # left ghost nodes
        u[m:m+3] = q_r  # right ghost nodes
        
        us[:, 0] = u[:m]
        
        # Integrate the PDE
        for i in range(1, snaps+1):
            for k in range(iters):
                for j in range(m):
                    # Compute C_hat to determine upwind direction
                    C_hat = 1 - u[j+1] - u[j-1]
                    
                    # Select upwind stencils based on sign of C_hat
                    if C_hat > 0:
                        # Positive characteristic velocity: upwind from left
                        ur = eno(u[j-2], u[j-1], u[j], u[j+1], u[j+2])
                        ul = eno(u[j-3], u[j-2], u[j-1], u[j], u[j+1])
                    else:
                        # Negative characteristic velocity: upwind from right
                        ur = eno(u[j+3], u[j+2], u[j+1], u[j], u[j-1])
                        ul = eno(u[j+2], u[j+1], u[j], u[j-1], u[j-2])
                    
                    # Compute final characteristic velocity
                    C_j = 1 - ur - ul
                    
                    # Update using traffic model
                    u1[j] = u[j] - (dt / dx) * C_j * (ur - ul)
                
                u[:m] = u1[:m]
                # Update ghost nodes
                u[-3:] = q_l
                u[m:m+3] = q_r
            
            us[:, i] = u[:m]
        
        # Plot the snapshots at t = 0, 0.5, 1.0, 1.5
        x = np.linspace(-1, 1, m)
        time_indices = [0, 50, 100, 150]
        times = [0, 0.5, 1.0, 1.5]
        
        plt.figure(figsize=(10, 6))
        for idx, t in zip(time_indices, times):
            plt.plot(x, us[:, idx], label=f't = {t}')
        
        plt.xlabel('x')
        plt.ylabel('Density')
        plt.title(case_name)
        plt.legend()
        plt.grid(True)
        plt.show()

    # Update parameters for traffic problem
    T = 1.5
    snaps = 150
    iters = int(T / dt) // snaps
    us = np.empty((m, snaps+1))  # Reallocate us with new snaps

    # Red light case
    print("Running red light case...")
    run_case(0.4, 1.0, "Red Light")

    # Green light case
    print("Running green light case...")
    run_case(1.0, 0.4, "Green Light")




def Problem2():
    # Define sinc function
    def sinc(x):
        if abs(x) < 1e-10:
            return 1.0
        return np.sin(np.pi * x) / (np.pi * x)

    # Define v(x) = sin(pi*x/2)
    def v(x):
        return np.sin(np.pi * x / 2)

    # Compute p_alpha(x)
    def p_alpha(x, alpha):
        result = 0.0
        for m in range(-alpha, alpha+1):
            x_m = m  # h = 1
            v_m = v(x_m)
            result += v_m * sinc(x - x_m)
        return result

    # Compute L2 error norm over [-5, 5]
    def compute_error(alpha):
        def integrand(x):
            diff = p_alpha(x, alpha) - v(x)
            return diff ** 2
        error_squared, _ = quad(integrand, -5, 5)
        return np.sqrt(error_squared)

    # Compute errors for alpha = 1, 2, 4, 8, ..., 512
    alphas = [2**k for k in range(10)]
    errors = []

    print("Computing errors for sinc interpolation...")
    with alive_bar(len(alphas)) as bar:
        for alpha in alphas:
            error = compute_error(alpha)
            errors.append(error)
            bar()

    # Fit to power law E_alpha = C * alpha^q using log-log regression
    log_alphas = np.log(alphas)
    log_errors = np.log(errors)
    coeffs = np.polyfit(log_alphas, log_errors, 1)
    q_fit = coeffs[0]
    C_fit = np.exp(coeffs[1])

    print(f"\nPower law fit: E_alpha = {C_fit:.6e} * alpha^{q_fit:.6f}")

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.loglog(alphas, errors, 'o-', label='Computed errors', markersize=8)
    plt.loglog(alphas, C_fit * np.array(alphas) ** q_fit, '--', label=f'Fit: $E_\\alpha = {C_fit:.3e} \\alpha^{{{q_fit:.3f}}}$')
    plt.xlabel('$\\alpha$')
    plt.ylabel('$E_\\alpha$')
    plt.title('L2 Error vs $\\alpha$ for Sinc Interpolation')
    plt.legend()
    plt.grid(True, which='both')
    plt.show()

    pass


def Problem3():
    # Part 2: Solve BVP for range of N

    # Exact solution u(x) = e^x(x^2-1)
    def u_exact(x):
        return np.exp(x) * (x**2 - 1)
    
    # Second derivative: u''(x) = e^x(x^2+4x+1)
    def u_double_prime(x):
        return np.exp(x) * (x**2 + 4*x + 1)
    
    # Source term f = u_xx + u^5
    def f_source(x):
        return u_double_prime(x) + u_exact(x)**5
    
    # Chebyshev nodes and differentiation matrices
    def chebyshev_nodes(N):
        return np.cos(np.pi * np.arange(N+1) / N)
    
    def chebyshev_diff_matrix(N):
        x = chebyshev_nodes(N)
        D = np.zeros((N+1, N+1))
        for i in range(N+1):
            for j in range(N+1):
                if i == j:
                    if i == 0:
                        D[i, j] = (2*N**2 + 1) / 6
                    elif i == N:
                        D[i, j] = -(2*N**2 + 1) / 6
                    else:
                        D[i, j] = -x[i] / (2 * (1 - x[i]**2))
                else:
                    c_i = 2 if (i == 0 or i == N) else 1
                    c_j = 2 if (j == 0 or j == N) else 1
                    D[i, j] = (c_i / c_j) * ((-1)**(i+j)) / (x[i] - x[j])
        return D
    
    # Solve BVP using Newton's method
    def solve_bvp(N):
        x = chebyshev_nodes(N)
        D = chebyshev_diff_matrix(N)
        D2 = D @ D
        f = np.array([f_source(xi) for xi in x])
        
        u = (x + 1)**2 * (x - 1)  # Initial guess
        
        for _ in range(50):
            residual = D2 @ u + u**5 - f
            residual[0] = u[0]
            residual[-1] = u[-1]
            
            if np.linalg.norm(residual) < 1e-10:
                break
            
            J = D2 + 5 * np.diag(u**4)
            J[0, :] = 0
            J[0, 0] = 1
            J[-1, :] = 0
            J[-1, -1] = 1
            
            u = u - np.linalg.solve(J, residual)
        
        return u, x
    
    # Compute errors for N = 4 to 64
    Ns = np.arange(4, 65, 2)
    errors = []
    
    for N in Ns:
        u_num, x = solve_bvp(N)
        u_ex = np.array([u_exact(xi) for xi in x])
        error = np.sqrt(np.sum((u_num - u_ex)**2) / (N+1))
        errors.append(error)
    
    # Plot semilog
    plt.figure(figsize=(10, 6))
    plt.semilogy(Ns, errors, 'o-', markersize=8)
    plt.xlabel('N')
    plt.ylabel('$\|p_N - u\|_2$')
    plt.title('L2 Error vs N for Chebyshev Spectral Method')
    plt.grid(True, which='both')
    plt.show()



if __name__ == "__main__":
    # Uncomment the problem you want to run
    # Problem1()
    # Problem2()
    # Problem3()
    pass