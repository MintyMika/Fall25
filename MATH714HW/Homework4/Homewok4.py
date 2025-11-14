import numpy as np
import matplotlib.pyplot as plt
from alive_progress import alive_bar

def HW4_Q2():
    # """
    # Dropping the last term in the Beam--Warming method from Eq.~\eqref{eq:beam_warming}
    #     gives
    #     \ begin{equation}
    #       U^{n+1}_j = U^n_j - \ frac{ak}{2h} (3 U^n_j - 4U^n_{j-1} + U^n_{j-2} ), \label{eq:sec_ord}
    #     \end{equation}
    #     which corresponds to forward Euler method in time, and a second-order
    #     one-sided derivative in space. Define $\ nu = ak/h$.
    #     \ begin{enumerate}[a)]
    #       \item Calculate the amplification factor $g(\xi)$ for a plane wave solution
    #             $U^0_j = e^{ijh\xi}$.
    #       \item Define $A(\xi)= |g(\xi)|^2$ and calculate a Taylor series for $A$ at
    #             $\xi=0$ up to second order. Using the Taylor series, explain why we
    #             consider the numerical scheme of Eq.~\eqref{eq:sec_ord} to be unstable
    #             regardless of the choice of timestep.
    #       \item Make two plots of $A(\xi)$ for $\nu = \ nicefrac1{100}$ using two different
    #             axis ranges:
    #             \ begin{itemize}
    #               \item $0\le h\xi \le 2\pi$ and $0.91 \le A \le 1.01$,
    #               \item $0\le h\xi \le 0.17$ and $1-10^{-6} \le A \le 1+10^{-6}$.
    #             \end{itemize}
    #       \item Write a program to simulate Eq.~\eqref{eq:sec_ord} on a periodic interval $[0,2\pi)$ using $N=40$ grid points
    #             and a grid spacing of $h = 2\pi/N$. Use the initial condition $u= \exp (2\sin x)$
    #             and $\ nu =\ nicefrac1{100}$. Plot the solution for $n=0,1000,2000,4000$. Define
    #             the root mean squared value of the solution,
    #             \ begin{equation}
    #               R(n) = \sqrt{\ frac{1}{N} \sum_{j=0}^{N-1} (U^n_j)^2}.
    #             \end{equation}
    #             Make a plot of $R$ over the range from $n=0$ to $n=10000$. You should find that $R$
    #             does not grow over time, indicating that the method is stable.
    #       \item Using the discrete Fourier transform, it can be shown that an arbitrary initial
    #             condition on the periodic interval can be written as
    #             \ begin{equation}
    #               U^0_j = \sum_{l=0}^{N-1} \ alpha_l e^{ijlh}
    #             \end{equation}
    #             for some constants $\ alpha_l$. Write down an expression for the general
    #             solution $U^n_j$. Using your answer, explain why your result in part (d)
    #             does not contradict the result in part (b).
    #     \end{enumerate}
    # """

    # The amplification factor is given by $g(\xi) = 1 - \frac{\nu}{2} (3 - 4e^{-ih\xi} + e^{-2ih\xi})$. 
    # The Taylor series expansion of $A(\xi) = |g(\xi)|^2$ about $\xi=0$ up to second order is $A(\xi) = 1 + \nu(1+\nu)(h\xi)^2 + O((h\xi)^4)$. Since the coefficient of the $(h\xi)^2$ term is positive for all $\nu > 0$, $A(\xi) > 1$ for small but nonzero $\xi$, indicating instability.

    # Now we can make the plots for part (c)
    nu = 1/100

    # Full range plot 0 <= hξ <= 2π
    theta = np.linspace(0, 2*np.pi, 4000)
    g = 1 - (nu/2) * (3 - 4*np.exp(-1j*theta) + np.exp(-2j*theta))
    A = np.abs(g)**2

    plt.figure(figsize=(7,4))
    plt.plot(theta, A, lw=1)
    plt.xlim(0, 2*np.pi)
    plt.ylim(0.91, 1.01)
    plt.xlabel('h xi')
    plt.ylabel('A(h xi)')
    plt.title('A(hxi) for ν = 1/100 (0 ≤ hxi ≤ 2 pi)')
    plt.grid(True)

    # Zoomed-in plot near 0: 0 <= hξ <= 0.17
    theta_zoom = np.linspace(0, 0.17, 2000)
    g_zoom = 1 - (nu/2) * (3 - 4*np.exp(-1j*theta_zoom) + np.exp(-2j*theta_zoom))
    A_zoom = np.abs(g_zoom)**2



    # Uncomment the following lines if you want to see the plots for (c) before part (d)
    # plt.figure(figsize=(7,4))
    # plt.plot(theta_zoom, A_zoom, lw=1)
    # plt.xlim(0, 0.17)
    # plt.ylim(1 - 1e-6, 1 + 1e-6)
    # plt.xlabel('h xi')
    # plt.ylabel('A(h xi)')
    # plt.title('Zoomed A(h xi) near 0 for ν = 1/100')
    # plt.grid(True)

    # plt.show() 


    # Now we can implement the numerical scheme for part (d)
    # Numerical simulation for part (d)
    N = 40
    h = 2 * np.pi / N
    x = np.linspace(0, 2*np.pi, N, endpoint=False)
    U = np.exp(2 * np.sin(x))

    max_n = 10000
    snap_times = {0, 1000, 2000, 4000}
    snapshots = {}
    R = np.zeros(max_n + 1)

    R[0] = np.sqrt(np.mean(U**2))
    if 0 in snap_times:
        snapshots[0] = U.copy()

    with alive_bar(max_n) as bar: # I totally thought this would take longer so I added a progress bar
        for n in range(1, max_n + 1):
            Um1 = np.roll(U, 1)   # U_{j-1}
            Um2 = np.roll(U, 2)   # U_{j-2}
            U = U - (nu/2.0) * (3*U - 4*Um1 + Um2)
            R[n] = np.sqrt(np.mean(U**2))
            if n in snap_times:
                snapshots[n] = U.copy()
            bar()

    # Plot snapshots at requested times
    plt.figure(figsize=(8,5))
    for t in sorted(snapshots):
        plt.plot(x, snapshots[t], label=f"n={t}")
    plt.xlabel("x")
    plt.ylabel("U")
    plt.title("Solution at n = 0, 1000, 2000, 4000")
    plt.legend()
    plt.grid(True)

    # Plot R(n) over time
    plt.figure(figsize=(8,4))
    plt.plot(np.arange(max_n + 1), R, lw=1)
    plt.xlabel("n")
    plt.ylabel("R(n)")
    plt.title("Root mean squared R(n) for n=0..10000")
    plt.grid(True)

    plt.show()
    pass




def HW4_Q3b():
    pass


def HW4_Q3c():
    pass

if __name__ == "__main__":
    HW4_Q2()

    pass