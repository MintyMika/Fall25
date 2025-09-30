def problem1():
    # Finite Difference Approximation Testing on the function f(x) = e^{-x}\tan(x). Also plot the log-log plot of error magnitude E vs H as well as a linear regression fit to E = C*H^p and determine C and p to three significant figures.
    # import statements
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
        

    # Define the function and its derivative
    def f(x):
        return np.exp(-x) * np.tan(x)

    def f_double_prime(x):
        return np.exp(-x) * (2 * np.tan(x)**3 - 2 * np.tan(x)**2 + 2 * np.tan(x) - 1)

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
        # \begin{equation*}
        # u(x) = -\frac{1}{2} e^{\pi} \cos x - \frac{1}{2} e^{\pi} \sin x - \frac{1}{2} e^{x}.
        # \end{equation*}
        return -(0.5 * np.exp(pi) * np.cos(x) + 0.5 * np.exp(pi) * np.sin(x) + 0.5 * np.exp(x))
    
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
        
    

    # # Optional  |Plotting the error trend
    # # Debugging |to make sure the solution looks correct
    # hs = [pi / n for n in ns]
    # plt.figure()
    # plt.loglog(hs, errors, marker='o')
    # plt.xlabel('h (log scale)')
    # plt.ylabel('E_n (log scale)')
    # plt.title('Error Trend (log-log)')
    # plt.grid(True, which="both", ls="--")
    # plt.show()

    # # Plot the equation and the numerical solution for n=160
    # U, x = solve_bvp(160)
    # u_ex = u_exact(x)
    # plt.figure()
    # plt.plot(x, U, label='Numerical Solution (n=160)')
    # plt.plot(x, u_ex, label='Exact Solution', linestyle='--')
    # plt.xlabel('x')
    # plt.ylabel('u(x)')
    # plt.title('Numerical vs Exact Solution')
    # plt.legend()
    # plt.grid()
    # plt.show()


    pass


def problem3():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

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
    
    def newtons_method(n, tol=10e-10, max_iter=100):
        h = 1 / n
        U = np.zeros(n - 1)  # Initial guess U^0 = 0
        for k in range(max_iter):
            # Evaluate the function F(U) at the current guess U
            F_val = F(U, h)
            
            # Evaluate the Jacobian matrix J_F(U) at the current guess U
            J_val = J_F(U, h)
            
            # Solve the linear system J_F(U) * delta_U = -F(U) to find the update delta_U
            delta_U = np.linalg.solve(J_val, -F_val)
            
            # Update the solution U with the computed delta_U
            U += delta_U
            
            # Check the termination condition: relative step size is less than the tolerance
            if np.linalg.norm(delta_U, 2) / np.linalg.norm(U, 2) < tol:
                print(f"Newton's method converged in {k+1} iterations.")
                break
        else:
            print("Newton's method did not converge within the maximum number of iterations.")
        
        # Return the computed solution U
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

    print(f"U[50] to three significant figures: {U_full[50]:.3}")


def problem4():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats


    def problem4_part_c():
        # Writing a program to solve the equation ∇²u = f on the triangular domain T with boundary conditions u(x) = 0 for x ∈ ∂T.

        s = np.sqrt(3) / 2
        def u_exact(x, y):
            return ((2 * y - np.sqrt(3))**2 - 3 * (2 * x - 1)**2) * np.sin(y)
        
        ns = [10, 20, 40, 80, 160]
        errors = []


        # From my derivation in part (a):
            #       \begin{equation*}
            #     \nabla^2_3 u_{i,j} = -\frac{4}{h^2} u_{i,j} + \frac{4}{3h^2} (u_{i+1,j} + u_{i,j-1} + u_{i-1,j+1} + u_{i-1,j} + u_{i,j+1} + u_{i+1,j-1}).
            #   \end{equation*}
        
        def finite_difference_approximation(u, i, j, h):
            return (-4 / h**2) * u[i, j] + (4 / (3 * h**2)) * (u[i + 1, j] + u[i, j - 1] + u[i - 1, j + 1] + u[i - 1, j] + u[i, j + 1] + u[i + 1, j - 1])

        def solve_system(n):
            h = 1 / n
            u = np.zeros((n, n))
            for j in range(1, n):
                for i in range(1, n - j):
                    u[i, j] = finite_difference_approximation(u, i, j, h)
            return u

        def error_measure(u_num, u_ex, n):
            h = 1 / n
            s = np.sqrt(3) / 2
            error_sum = 0
            for j in range(1, n):
                for i in range(1, n - j):
                    error_sum += (u_num[i, j] - u_ex(i * h + 0.5 * j * h, j * h))**2
            return np.sqrt((s / (2 * n**2)) * error_sum)

        for n in ns:
            u_num = np.zeros((n, n))
            h = 1 / n
            u_num = solve_system(n)
            errors.append(error_measure(u_num, u_exact, n))



        hs = [1 / n for n in ns]
        plt.figure()
        plt.loglog(hs, errors, marker='o')
        plt.xlabel('h')
        plt.ylabel('Error')
        plt.title('Error Convergence')
        plt.grid()
        plt.show()

        # Linear regression
        log_h = np.log(hs)
        log_errors = np.log(errors)
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_h, log_errors)
        C = np.exp(intercept)
        p = slope

        # Print C and p
        print(f"C: {C:.3}, p: {p:.3}")


    def problem4_part_d():
        # Writing a program to solve the equation ∇²u = f on the triangular domain T with boundary conditions u(x) = 0 for x ∈ ∂T.

        s = np.sqrt(3) / 2
        def u_exact(x, y):
            return ((2 * y - np.sqrt(3))**2 - 3 * (2 * x - 1)**2) * np.sin(y)
        
        ns = [10, 20, 40, 80, 160]
        errors = []


        # From my derivation in part (a):
            #       \[
            #     c_0 = -\frac{4}{h^2}, \quad c_1 = c_2 = c_3 = c_4 = c_5 = c_6 = \frac{2}{3h^2}.
            #   \]
        
        def finite_difference_approximation(u, i, j, h):
            return (-4 / h**2) * u[i, j] + (2 / (3 * h**2)) * (u[i + 1, j] + u[i, j - 1] + u[i - 1, j + 1] + u[i - 1, j] + u[i, j + 1] + u[i + 1, j - 1])

        def solve_system(n):
            h = 1 / n
            u = np.zeros((n, n))
            for j in range(1, n):
                for i in range(1, n - j):
                    u[i, j] = finite_difference_approximation(u, i, j, h)
            return u

        def error_measure(u_num, u_ex, n):
            h = 1 / n
            s = np.sqrt(3) / 2
            error_sum = 0
            for j in range(1, n):
                for i in range(1, n - j):
                    error_sum += (u_num[i, j] - u_ex(i * h + 0.5 * j * h, j * h))**2
            return np.sqrt((s / (2 * n**2)) * error_sum)

        for n in ns:
            u_num = np.zeros((n, n))
            h = 1 / n
            u_num = solve_system(n)
            errors.append(error_measure(u_num, u_exact, n))



        hs = [1 / n for n in ns]
        plt.figure()
        plt.loglog(hs, errors, marker='o')
        plt.xlabel('h')
        plt.ylabel('Error')
        plt.title('Error Convergence')
        plt.grid()
        plt.show()

        # Linear regression
        log_h = np.log(hs)
        log_errors = np.log(errors)
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_h, log_errors)
        C = np.exp(intercept)
        p = slope

        # Print C and p
        print(f"C: {C:.3}, p: {p:.3}")


    problem4_part_c()
    problem4_part_d()

    pass

def main():
    print("Uncomment the problem you want to run in the main function.")
    # problem1()
    # problem2()
    # problem3() 
    # problem4()

if __name__ == "__main__":
    main()