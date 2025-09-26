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
    # Problem Statement 2:




    
    pass

def problem3():
    # TODO: Implement solution for Homework Problem 3
    pass

def problem4():
    # TODO: Implement solution for Homework Problem 4
    pass

def main():
    problem1()
    problem2()
    problem3()
    problem4()

if __name__ == "__main__":
    main()