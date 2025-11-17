import numpy as np
import matplotlib.pyplot as plt
from alive_progress import alive_bar

def HW4_Q2():
    # Part (c)
    nu = 1/100

    # Full range plot 0 <= h xi <= 2 pi
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

    # Zoomed-in plot near 0: 0 <= h xi <= 0.17
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




def HW4_Q3():
    # We first note that the maximum value of A(x) occurs at x = pi/2, giving A_max = 2 + 4/3 = 10/3. Thus, c = 10/3.

    # Now we can implement the Lax-Wendroff scheme
    # Define A(x) and CFL speed
    def A_of(x):
        return 2.0 + (4.0/3.0) * np.sin(x)

    c = np.max(np.abs(A_of(np.linspace(0, 2*np.pi, 1000))))

    # Parameters for part (b)
    T = 3.0 * np.pi / np.sqrt(5.0)

    def run_lax_wendroff(m, q0, dt_factor=1/3):
        h = 2*np.pi / m
        x_centers = np.linspace(h/2, 2*np.pi - h/2, m)
        A_centers = A_of(x_centers)             # A_i = A((i+1/2)h)
        x_edges = np.linspace(0, 2*np.pi, m, endpoint=False)
        A_edges = A_of(x_edges)                 # A_{i-1/2} = A(i h)

        dt = h / (dt_factor * c)                # note dt_factor=3 => dt = h/(3c) when dt_factor=3
        # To hit t=T exactly, adjust number of steps
        nsteps = int(round(T / dt))
        dt = T / nsteps

        Q = q0.copy()
        snapshot_indices = [0, int(round(nsteps/4)), int(round(nsteps/2)), int(round(3*nsteps/4)), nsteps]
        snapshots = {}

        for n in range(nsteps+1):
            if n in snapshot_indices:
                snapshots[n] = Q.copy()
            if n == nsteps:
                break

            Q_im1 = np.roll(Q, 1)
            A_im1 = np.roll(A_centers, 1)

            F_imhalf = 0.5 * (A_im1 * Q_im1 + A_centers * Q) - (A_edges * dt) / (2*h) * (A_centers * Q - A_im1 * Q_im1)
            F_iphalf = np.roll(F_imhalf, -1)

            Q = Q - (dt / h) * (F_iphalf - F_imhalf)

        return x_centers, snapshots, Q

    # Part (b): m = 512
    m = 512
    h = 2*np.pi / m
    x_centers = np.linspace(h/2, 2*np.pi - h/2, m)
    q0_smooth = np.exp(np.sin(x_centers) + 0.5 * np.sin(4*x_centers))

    # Note: dt_factor used here as 3 to get dt = h/(3c)
    _, snaps_smooth, Q_T_smooth = run_lax_wendroff(m, q0_smooth, dt_factor=3)

    plt.figure(figsize=(8,5))
    for n in sorted(snaps_smooth):
        plt.plot(x_centers, snaps_smooth[n], label=f"n={n}")
    plt.xlabel("x")
    plt.ylabel("q")
    plt.title("Lax--Wendroff snapshots, m=512")
    plt.legend()
    plt.grid(True)

    # Part (c)
    m_list = [256, 512, 1024, 2048]
    errors = []
    h_list = []
    for mtest in m_list:
        h_test = 2*np.pi / mtest
        x_ct, _, Q_T = run_lax_wendroff(mtest, np.exp(np.sin(np.linspace(h_test/2, 2*np.pi-h_test/2, mtest)) + 0.5*np.sin(4*np.linspace(h_test/2, 2*np.pi-h_test/2, mtest))), dt_factor=3)
        q0_test = np.exp(np.sin(x_ct) + 0.5*np.sin(4*x_ct))   # exact at T equals initial
        err = np.sqrt(h_test * np.sum((Q_T - q0_test)**2))
        errors.append(err)
        h_list.append(h_test)

    # Estimate convergence order
    rate = np.polyfit(np.log(h_list), np.log(errors), 1)[0]

    plt.figure(figsize=(6,4))
    plt.loglog(h_list, errors, '-o')
    plt.xlabel('h')
    plt.ylabel('L2 error at T')
    plt.title(f'Convergence, estimated order ≈ {abs(rate):.2f}')
    plt.grid(True, which="both")

    # Part (d)
    def alternate_ic(x):
        return np.maximum(np.pi/2 - np.abs(x - np.pi), 0.0)

    # m=512 snapshots
    m = 512
    h = 2*np.pi / m
    x_centers = np.linspace(h/2, 2*np.pi - h/2, m)
    q0_alternate = alternate_ic(x_centers)
    _, snaps_alternate, Q_T_alternate = run_lax_wendroff(m, q0_alternate, dt_factor=3)

    plt.figure(figsize=(8,5))
    for n in sorted(snaps_alternate):
        plt.plot(x_centers, snaps_alternate[n], label=f"n={n}")
    plt.xlabel("x")
    plt.ylabel("q")
    plt.title("Lax--Wendroff snapshots Alternative initial condition, m=512")
    plt.legend()
    plt.grid(True)

    # Convergence
    errors_alternate = []
    h_list_alternate = []
    for mtest in m_list:
        h_test = 2*np.pi / mtest
        x_ct = np.linspace(h_test/2, 2*np.pi - h_test/2, mtest)
        q0_test = alternate_ic(x_ct)
        _, _, Q_T = run_lax_wendroff(mtest, q0_test, dt_factor=3)
        err = np.sqrt(h_test * np.sum((Q_T - q0_test)**2))
        errors_alternate.append(err)
        h_list_alternate.append(h_test)

    rate_alternate = np.polyfit(np.log(h_list_alternate), np.log(errors_alternate), 1)[0]

    plt.figure(figsize=(6,4))
    plt.loglog(h_list_alternate, errors_alternate, '-o')
    plt.xlabel('h')
    plt.ylabel('L2 error at T')
    plt.title(f'Convergence Alternative initial condition, estimated order ≈ {abs(rate_alternate):.3f}')
    plt.grid(True, which="both")

    plt.show()

    pass


if __name__ == "__main__":
    # HW4_Q2() # Uncomment to run Q2
    HW4_Q3() # Uncomment to run Q3
    pass