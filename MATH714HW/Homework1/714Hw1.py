def problem1():
    # Problem statement: 
    # Let $f: \R \to \R$ be a smooth function, and let $x_1<x_2<x_3<x_4$ be
    #     four increasing values.
    #     \begin{enumerate}[(a)]
    #       \item Derive a finite difference (FD) approximation for $f''(x_2)$ that
    #             is as accurate as possible, based on the four values of $f_1=f(x_1),
    #               \ldots, f_4 = f(x_4)$. Calculate an expression for the dominant term
    #             in the error.
    #       \item Write a program to test the FD approximation
    #             on the function
    #             \begin{equation}
    #               f(x) = e^{-x} \tan x.
    #             \end{equation}
    #             Consider step sizes of $H=10^{-k/100}$ for $k={100,101,\ldots,300}$.
    #             For each $H$, set $x_1=0$ and $x_4=H$. Choose $x_2$ and $x_3$ as uniformly
    #             randomly distributed random numbers over the range from $0$ to
    #             $H$.\footnote{If $x_2>x_3$, then swap the two values to ensure the
    #               ordering is preserved. If $x_2=x_3$, then choose new random numbers.}
    #             Make a log--log plot showing the absolute error magnitude $E$ of the FD
    #             approximation versus $H$. Use linear regression to fit the data to
    #             \begin{equation}
    #               E = C H^p
    #             \end{equation}
    #             and determine $C$ and $p$ to three significant figures.\footnote{Since sample
    #               points for the FD approximation are randomly chosen, there will be
    #               small variations in the values of $C$ and $p$ that you compute.}

    # Finite Difference Approximation Testing on the function f(x) = e^{-x}\tan(x). Also plot the log-log plot of error magnitude E vs H as well as a linear regression fit to E = C*H^p and determine C and p to three significant figures.
    # import statements
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from math import tan, exp, pi

    # Define the function and its derivative
    def f(x):
        return exp(-x) * tan(x)


    pass

def problem2():
    # TODO: Implement solution for Homework Problem 2
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