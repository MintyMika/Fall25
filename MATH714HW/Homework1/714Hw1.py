def problem1():
    # Finite Difference Approximation Testing on the function f(x) = e^{-x}\tan(x). Also plot the log-log plot of error magnitude E vs H as well as a linear regression fit to E = C*H^p and determine C and p to three significant figures.
    # import statements
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from math import tan, exp, pi

    # Define the function and its derivative
    def f(x):
        return exp(-x) * tan(x)

    def f_double_prime(x):
        return exp(-x) * (2 * tan(x)**3 - 2 * tan(x)**2 + 2 * tan(x) - 1)
    
    # Finite difference approximation for f''(x2)
    def finite_difference_approximation(f1, f2, f3, f4, x1, x2, x3, x4):
        # Using the method of undetermined coefficients I derived the following formula:
        a = (2*x2 - x3 - x4) / ((x1 - x2)*(x1 - x3)*(x1 - x4))
        b = (2*x2 - x1 - x4) / ((x2 - x1)*(x2 - x3)*(x2 - x4))
        c = (2*x2 - x1 - x4) / ((x3 - x1)*(x3 - x2)*(x3 - x4))
        d = (2*x2 - x1 - x3) / ((x4 - x1)*(x4 - x2)*(x4 - x3))

        return a*f1 + b*f2 + c*f3 + d*f4
    
    # Parameters
    H_values = [10**(-k/100) for k in range(100, 301)]  #
    errors = []
    x1 = 0
    x2 = None  # Will be randomly chosen < x3
    x3 = None  # Will be randomly chosen > x2
    x4 = None  # Will be set to H in the loop
    true_value = f_double_prime(0)  # f''(0)

    for H in H_values:
        x4 = H
        # Randomly choose x2 and x3
        while True:
            x2, x3 = np.sort(np.random.uniform(0, H, 2))
            if x2 != x3:
                break
        
        # Make sure x2 < x3
        x2, x3 = min(x2, x3), max(x2, x3)


        # Function values
        f1 = f(x1)
        f2 = f(x2)
        f3 = f(x3)
        f4 = f(x4)

        approx = finite_difference_approximation(f1, f2, f3, f4, x1, x2, x3, x4)
        error = abs(approx - true_value)
        errors.append(error)

    # Log-log plot
    plt.figure()
    plt.loglog(H_values, errors, marker='o', linestyle='None')
    plt.xlabel('H')
    plt.ylabel('Error')
    plt.title('Error of Finite Difference Approximation')
    plt.grid()

    # Linear regression
    log_H = np.log(H_values)
    log_errors = np.log(errors)
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_H, log_errors)
    C = np.exp(intercept)
    p = slope

    # Print C and p
    print(f"C: {C:.3}, p: {p:.3}")

    plt.loglog(H_values, C * H_values**p, linestyle='--', color='r')
    plt.legend(['Error', f'Fit: E = {C:.3f} H^{p:.3f}'])
    plt.show()


def problem2():

    # TODO: Fix the matrix A as the boundary conditions are not correctly implemented.

    # Problem Statement 2:
    # \textbf{Mixed boundary value problem (5 points).}
    #     For a smooth function $u(x)$ and source term $f(x)$, consider the two
    #     point boundary value problem (BVP)
    #     \begin{equation}
    #       u''+u = f(x)
    #     \end{equation}
    #     on the domain $x\in [0,\pi]$, using the mixed boundary conditions
    #     \begin{equation}
    #       u'(0) - u(0) = 0, \qquad u'(\pi) + u(\pi) =0.
    #     \end{equation}
    #     \begin{enumerate}[(a)]
    #       \item Use a mesh width of $h=\pi/n$ where $n\in \N$ and introduce grid
    #             points at $x_j=jh$ for $j=0,1,\ldots, n$. Construct a second-order
    #             accurate finite-difference method for this BVP. Write
    #             your method as a linear system of the form $AU=F$.
    #       \item Construct the exact solution $u(x)$ to the BVP when $f(x)=-e^x$.
    #       \item Verify that your method is second-order accurate by solving the BVP
    #             with $f(x)=-e^x$ using $n=20,40,80,160$. For each $n$, construct the
    #             error measure
    #             \begin{equation}
    #               E_n = \sqrt{ h \sum_{j=0}^n q_j (U_j - u(x_j))^2 }
    #             \end{equation}
    #             where $U_j$ is the numerical solution at $x_j$. Here, $q_j = \tfrac12$
    #             when $j\in \{0,n\}$ and $q_j=1$ otherwise. Present your results in a
    #             table, and comment on whether the trend in the errors is expected for a
    #             second-order method.
    #     \end{enumerate}



    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from math import tan, exp, pi

    def f(x):
        return -np.exp(x)
    
    def error_measure(U, u_exact, h):
        n = len(U) - 1
        q = np.ones(n + 1)
        q[0] = q[n] = 0.5
        return np.sqrt(h * np.sum(q * (U - u_exact)**2))
    
    def u_exact(x):
        return -0.5 * np.exp(x) + 0.5 * (np.cos(x) + np.sin(x)) / np.cos(pi)
    
    def solve_bvp(n):
        h = pi / n
        x = np.linspace(0, pi, n + 1)
        A = np.zeros((n + 1, n + 1))
        F = f(x)

        # Interior points
        for j in range(1, n):
            A[j, j - 1] = 1 / h**2
            A[j, j] = -2 / h**2 + 1
            A[j, j + 1] = 1 / h**2

        # Boundary conditions
        # At x=0: Use forward difference for u'(0): u'(0) ≈ (-3U[0] + 4U[1] - U[2])/(2h)
        # u'(0) - u(0) = 0 => [(-3U[0] + 4U[1] - U[2])/(2h)] - U[0] = 0
        # => (-3U[0] + 4U[1] - U[2]) - 2h*U[0] = 0
        # => (-3 - 2h)U[0] + 4U[1] - U[2] = 0
        A[0, 0] = (-3 - 2*h) #/ (2*h)
        A[0, 1] = 4 #/ (2*h)
        A[0, 2] = -1 #/ (2*h)
        F[0] = 0

        # At x=pi: Use backward difference for u'(pi): u'(pi) ≈ (3U[n] - 4U[n-1] + U[n-2])/(2h)
        # u'(pi) + u(pi) = 0 => [(3U[n] - 4U[n-1] + U[n-2])/(2h)] + U[n] = 0
        # => (3U[n] - 4U[n-1] + U[n-2]) + 2h*U[n] = 0
        # => (3 + 2h)U[n] - 4U[n-1] + U[n-2] = 0
        A[n, n] = (3 + 2*h) #/ (2*h)
        A[n, n - 1] = -4 #/ (2*h)
        A[n, n - 2] = 1 #/ (2*h)
        F[n] = 0

        # Put the above back in later. Use the regular second-order finite difference for now.
        # A[0, 0] = (-(2*(1 + h)) / h**2) + 1
        # A[0, 1] = 2 / h**2
        # F[0] = f(0)
        # A[n, n] = (-(2*(1 + h)) / h**2) + 1
        # A[n, n - 1] = 2 / h**2
        # F[n] = f(pi) 

        # TODO: Fix this as the matrix A is incorrect because the loglog plot is not linear




        # Solve the linear system
        U = np.linalg.solve(A, F)
        return U, x

    ns = [20, 40, 80, 160]
    errors = []
    print(f"{'n':<10}{'E_n':<20}")
    for n in ns:
        U, x = solve_bvp(n)
        u_ex = u_exact(x)
        h = pi / n
        E_n = error_measure(U, u_ex, h)
        errors.append(E_n)
        print(f"{n:<10}{E_n:<20.10e}")
        # The trend in the errors should show that as n increases (and thus h decreases), the error E_n decreases at a rate consistent with second-order accuracy, i.e., E_n should be proportional to h^2.
        # This can be verified by observing the ratio of errors for successive n values.

    # Optional: Plotting the error trend
    # hs = [pi / n for n in ns]
    plt.figure()
    plt.loglog(ns, errors, marker='o')
    plt.xlabel('n (log scale)')
    plt.ylabel('E_n (log scale)')
    plt.title('Error Trend (log-log)')
    plt.grid(True, which="both", ls="--")
    plt.show()

    # Plot the equation and the numerical solution for n=160
    U, x = solve_bvp(160)
    u_ex = u_exact(x)
    plt.figure()
    plt.plot(x, U, label='Numerical Solution (n=160)')
    plt.plot(x, u_ex, label='Exact Solution', linestyle='--')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Numerical vs Exact Solution')
    plt.legend()
    plt.grid()
    plt.show()

    # Based on the regular plot of the numerical and exact solution, the numerical solution should closely follow the exact solution, indicating the accuracy of the finite difference method but it doesn't look like the function f(x) = -e^x 


    # Present your results in a table, and comment on whether the trend in the errors is expected for a second-order method.
    # The graph looks like the function f(x) = \sqrt{x} is it supposed to look like that? 

    pass


def problem3():
    # Problem Statement 3:
    # Consider the nonlinear BVP
    #     \begin{equation}
    #       u''(x) - 80\cos u(x) = 0
    #     \end{equation}
    #     with the boundary conditions $u(0)=0$ and $u(1)=10$. Define the mesh width
    #     $h=\tfrac1n$ for $n\in \N$, and introduce gridpoints $x_j=jh$ for
    #     $j=0,1,\ldots,n$.

    #     Let $U_i$ be the approximation of $u(x_i)$. From the boundary conditions,
    #     $U_0=0$ and $U_n=10$. Let $U=(U_1,U_2,\ldots,U_{n-1})$ be the vector of
    #     unknown function values, and write $F(U)=0$ as the nonlinear system
    #     of algebraic equations from the finite-difference approximation, where
    #     $F:\R^{n-1} \to \R^{n-1}$. The components are
    #     \begin{align}
    #       F_1(U)     & = \frac{U_2 - 2U_1}{h^2} - 80 \cos U_1,                                                 \\
    #       F_i(U)     & = \frac{U_{i+1} - 2U_i+ U_{i-1}}{h^2} - 80 \cos U_i \qquad \text{for $i=2,\ldots,n-2$,} \\
    #       F_{n-1}(U) & = \frac{10 - 2U_{n-1} + U_{n-2}}{h^2} - 80 \cos U_{n-1}.
    #     \end{align}
    #     \begin{enumerate}[a)]
    #       \item Calculate the Jacobian $J_F \in \R^{(n-1)\times(n-1)}$ for the function $F$ and
    #             describe its structure.
    #       \item Use Newton's method to solve the BVP for $n=100$, using the Jacobian matrix from part (a). Write $U^k$ to indicate the $k$th Newton step, and start with the initial guess $U^0=0$. Terminate Newton's method when the relative step size $\|\Delta U^k\|_2 / \|U^k\|_2$ is less than $10^{-10}$. Plot the solution $U$ over the interval $[0,1]$, and report the value of $U_{50}$ to three significant figures.
    #     \end{enumerate}
    # the Jacobian J_F is a tridiagonal matrix 

    # the Jacobian matrix $J_F$ has the following tridiagonal structure:
    #     \begin{equation*}
    #       J_F = \begin{pmatrix}
    #         -\frac{2}{h^2} + 80 \sin U_1 & \frac{1}{h^2} & 0 & \cdots & 0 \\
    #         \frac{1}{h^2} & -\frac{2}{h^2} + 80 \sin U_2 & \frac{1}{h^2} & \cdots & 0 \\
    #         0 & \frac{1}{h^2} & -\frac{2}{h^2} + 80 \sin U_3 & \cdots & 0 \\
    #         \vdots & \vdots & \vdots & \ddots & \vdots \\
    #         0 & 0 & 0 & \cdots & -\frac{2}{h^2} + 80 \sin U_{n-1}
    #       \end{pmatrix}.
    #     \end{equation*}

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from math import tan, exp, pi

    def F(U, h):
        n = len(U) + 1
        F_vec = np.zeros(n - 1)
        F_vec[0] = (U[1] - 2 * U[0]) / h**2 - 80 * np.cos(U[0])
        for i in range(1, n - 2):
            F_vec[i] = (U[i + 1] - 2 * U[i] + U[i - 1]) / h**2 - 80 * np.cos(U[i])
        F_vec[n - 2] = (10 - 2 * U[n - 2] + U[n - 3]) / h**2 - 80 * np.cos(U[n - 2])
        return F_vec
    
    def J_F(U, h):
        n = len(U) + 1
        J = np.zeros((n - 1, n - 1))
        J[0, 0] = -2 / h**2 + 80 * np.sin(U[0])
        J[0, 1] = 1 / h**2
        for i in range(1, n - 2):
            J[i, i - 1] = 1 / h**2
            J[i, i] = -2 / h**2 + 80 * np.sin(U[i])
            J[i, i + 1] = 1 / h**2
        J[n - 2, n - 3] = 1 / h**2
        J[n - 2, n - 2] = -2 / h**2 + 80 * np.sin(U[n - 2])
        return J
    
    def newtons_method(n, tol=1e-10, max_iter=100):
        h = 1 / n
        U = np.zeros(n - 1)  # Initial guess U^0 = 0
        for k in range(max_iter):
            F_val = F(U, h)
            J_val = J_F(U, h)
            delta_U = np.linalg.solve(J_val, -F_val)
            U += delta_U
            if np.linalg.norm(delta_U, 2) / np.linalg.norm(U, 2) < tol:
                break
        return U
    
    n = 100
    U = newtons_method(n)
    x = np.linspace(0, 1, n + 1)
    U_full = np.zeros(n + 1)
    U_full[1:n] = U
    U_full[n] = 10  # Boundary condition at x=1

    # Plot the solution
    plt.figure()
    plt.plot(x, U_full, label='Numerical Solution')
    plt.xlabel('x')
    plt.ylabel('U')
    plt.title('Solution of BVP')
    plt.legend()
    plt.grid()
    plt.show()

    print(f"U[50] to three significant figures: {U_full[50]:.3f}")

    pass

def problem4():
    # TODO: Implement solution for Homework Problem 4
    pass

def main():
    # problem1()
    # problem2() # TODO: Fix the matrix A as the boundary conditions are not correctly implemented.
    # problem3() # TODO: Look over this and veryify it is correct.
    # problem4()

if __name__ == "__main__":
    main()