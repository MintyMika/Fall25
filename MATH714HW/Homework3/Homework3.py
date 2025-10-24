import numpy as np 
import matplotlib.pyplot as plt
from alive_progress import alive_bar

def question_1():
    # We have calculated that the stability function the given RK method is R(z) = 1 + z + z^2/2 and that using the boundary locus method the function becomes z(\theta) = -1 \pm \sqrt{2e^{i\theta} - 1}
    # To plot the stability region, we can parametrize theta from 0 to 2pi and compute z(theta)
    theta = np.linspace(0, 2*np.pi, 400)
    z1 = -1 + np.sqrt(2*np.exp(1j*theta) - 1)
    z2 = -1 - np.sqrt(2*np.exp(1j*theta) - 1)

    plt.figure(figsize=(6,6))
    plt.plot(z1.real, z1.imag, 'b-', label='Branch 1')
    plt.plot(z2.real, z2.imag, 'r-', label='Branch 2')
    plt.fill_between(np.concatenate([z1.real, z2.real[::-1]]),
                     np.concatenate([z1.imag, z2.imag[::-1]]),
                     color='lightgray', alpha=0.5)
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(0, color='black', lw=0.5)
    plt.title('Region of Absolute Stability for the RK Method')
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    plt.show()

    # Code below can be uncommented to plot |R(z)| <= 1 instead of using the boundary locus method

    # # Alternatively, we can plot |R(z)| <= 1
    # def R(z):
    #     return 1 + z + z**2 / 2
    
    # # Plotting on a -4 to 2 real and -3 to 3 imaginary grid
    # x = np.linspace(-4, 2, 400)
    # y = np.linspace(-3, 3, 400)
    # X, Y = np.meshgrid(x, y)
    # Z = X + 1j*Y
    # RZ = R(Z)
    # stability_region = np.abs(RZ) <= 1
    # plt.figure(figsize=(6,6))
    # plt.contourf(X, Y, stability_region, levels=[0, 0.5, 1], colors=['lightgray', 'lightblue'])
    # plt.axhline(0, color='black', lw=0.5)
    # plt.axvline(0, color='black', lw=0.5)
    # plt.title('Region of Absolute Stability for the RK Method')
    # plt.xlabel('Re(z)')
    # plt.ylabel('Im(z)')
    # plt.grid(True)
    # plt.axis('equal')
    # plt.show()

    pass

def question_2():
    # I decided not to use the boundary locus method and instead directly plot |R(z)| <= 1
    # Plot |R(z)| <= 1
    def R(z):
        return (5*z + 12) / ((z - 4)*(z - 3))
    
    # Plotting on a -20 to 20 real and -20 to 20 imaginary grid

    x = np.linspace(-20, 20, 400)
    y = np.linspace(-20, 20, 400)
    # Build a 2D grid of (x, y) coordinates over the complex plane
    X, Y = np.meshgrid(x, y)
    # Form the complex grid Z = X + iY (each point is a complex number)
    Z = X + 1j*Y
    # Evaluate the stability function R at every point on the complex grid
    RZ = R(Z)

    stability_region = np.abs(RZ) <= 1
    plt.figure(figsize=(6,6))
    plt.contourf(X, Y, stability_region, levels=[0, 0.5, 1], colors=['lightgray', 'lightblue'])
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(0, color='black', lw=0.5)
    plt.title('Region of Absolute Stability for the TR-BDF2 Method')
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')
    plt.grid(True)
    plt.axis('equal')
    # Adding a legend for the filled contour
    plt.plot([], [], 'lightblue', label='|R(z)| ≤ 1')
    plt.plot([], [], 'lightgray', label='|R(z)| > 1')
    plt.legend()

    plt.show()

    pass


def question_4():
    # Load the map data
    path_to_file = r"MATH714HW/Homework3/714_hw3_files/van_vleck.txt"
    # Load the map
    map_data = np.loadtxt(path_to_file)
    
    # Problem parameters
    h = 0.225  # grid spacing in meters (22.5 cm)
    b = 0.55   # diffusion coefficient in m^2/s
    k = h**2 / (6 * b)  # timestep
    
    # Grid dimensions
    n_rows, n_cols = map_data.shape  # 73 x 160

    # For verification of map dimensions
    # print(f"Map dimensions: {n_rows} rows x {n_cols} cols")
    
    # Initialize concentration field
    u = np.zeros((n_rows, n_cols))
    
    # Pizza region P: 36 <= i < 40, 44 <= j < 48
    pizza_i = slice(36, 40)
    pizza_j = slice(44, 48)
    u[pizza_i, pizza_j] = 1.0
    
    # Professor locations (i, j)
    prof_T = (31, 14)
    prof_Q = (58, 103)
    prof_C = (58, 147)
    
    # Track when professors are distracted
    distraction_times = {'T': None, 'Q': None, 'C': None}
    distraction_threshold = 1e-4
    
    # Track concentration at professor locations over time
    times = []
    conc_T = []
    conc_Q = []
    conc_C = []
    
    # Target times for plotting
    plot_times = [1.0, 5.0, 25.0, 100.0]
    plot_fields = {}
    
    # Time stepping
    t = 0.0
    max_time = 110.0
    
    print(f"Timestep k = {k:.6f} s")
    print(f"Starting simulation to t = {max_time} s...")
    
    step = 0
    while t <= max_time:
        # Store data for plotting
        if step % 100 == 0:
            times.append(t)
            conc_T.append(u[prof_T])
            conc_Q.append(u[prof_Q])
            conc_C.append(u[prof_C])
        
        # Check for distraction times
        if distraction_times['T'] is None and u[prof_T] > distraction_threshold:
            distraction_times['T'] = t
            print(f"Professor T distracted at t = {t:.1f} s")
        if distraction_times['Q'] is None and u[prof_Q] > distraction_threshold:
            distraction_times['Q'] = t
            print(f"Professor Q distracted at t = {t:.1f} s")
        if distraction_times['C'] is None and u[prof_C] > distraction_threshold:
            distraction_times['C'] = t
            print(f"Professor C distracted at t = {t:.1f} s")
        
        # Store fields at target times
        for target_t in plot_times:
            if target_t not in plot_fields and t >= target_t:
                plot_fields[target_t] = u.copy()
                print(f"Saved field at t = {target_t} s")
        
        # Update concentration field using finite difference
        u_new = u.copy()
        
        for i in range(n_rows):
            for j in range(n_cols):
                # Skip if it's a wall
                if map_data[i, j] == 1:
                    continue
                
                # Skip if it's in the pizza region
                if 36 <= i < 40 and 44 <= j < 48:
                    continue
                
                # Get neighbors, applying ghost node boundary conditions
                u_im1 = u[i-1, j] if i > 0 and map_data[i-1, j] == 0 else u[i, j]
                u_ip1 = u[i+1, j] if i < n_rows-1 and map_data[i+1, j] == 0 else u[i, j]
                u_jm1 = u[i, j-1] if j > 0 and map_data[i, j-1] == 0 else u[i, j]
                u_jp1 = u[i, j+1] if j < n_cols-1 and map_data[i, j+1] == 0 else u[i, j]
                
                # Apply finite difference scheme
                laplacian = (u_ip1 + u_jm1 - 4*u[i, j] + u_im1 + u_jp1) / (h**2)
                u_new[i, j] = u[i, j] + k * b * laplacian
        
        u = u_new
        
        # Keep pizza region fixed
        u[pizza_i, pizza_j] = 1.0
        
        t += k
        step += 1
        
        # Printed progress every 1000 steps to make sure it's running
        # if step % 1000 == 0:
        #     print(f"Step {step}, t = {t:.2} s")
    
    # Part (c): Print distraction times
    print("\n=== Part (c) Results ===")
    print(f"Professor T distracted at: {distraction_times['T']:.1f} s")
    print(f"Professor Q distracted at: {distraction_times['Q']:.1f} s")
    print(f"Professor C distracted at: {distraction_times['C']:.1f} s")
    
    # Part (b): Plot concentration fields at target times
    for target_t in sorted(plot_fields.keys()):
        field = plot_fields[target_t]
        plt.figure(figsize=(10, 5))
        plt.imshow(field**0.25, cmap='hot', origin='upper', extent=[0, n_cols, n_rows, 0])
        plt.imshow(map_data, cmap='gray', alpha=0.3, origin='upper', extent=[0, n_cols, n_rows, 0])
        plt.colorbar(label='[u]^(1/4)')
        plt.title(f'Scaled Smell Concentration at t = {target_t} s')
        plt.xlabel('j (grid points)')
        plt.ylabel('i (grid points)')
        plt.plot([pizza_j.start, pizza_j.stop-1, pizza_j.stop-1, pizza_j.start, pizza_j.start],
                [pizza_i.start, pizza_i.start, pizza_i.stop-1, pizza_i.stop-1, pizza_i.start],
                'b-', linewidth=2, label='Pizza')
        plt.plot(prof_T[1], prof_T[0], 'go', markersize=10, label='Prof T')
        plt.plot(prof_Q[1], prof_Q[0], 'mo', markersize=10, label='Prof Q')
        plt.plot(prof_C[1], prof_C[0], 'co', markersize=10, label='Prof C')
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    # Part (d): Semilog plot of concentration at professor locations
    plt.figure(figsize=(10, 6))
    plt.semilogy(times, np.maximum(conc_T, 1e-10), 'g-', label='Professor T', linewidth=2)
    plt.semilogy(times, np.maximum(conc_Q, 1e-10), 'm-', label='Professor Q', linewidth=2)
    plt.semilogy(times, np.maximum(conc_C, 1e-10), 'c-', label='Professor C', linewidth=2)
    plt.axhline(distraction_threshold, color='r', linestyle='--', label=f'Distraction threshold ({distraction_threshold})')
    plt.ylim(1e-10, 1)
    plt.xlim(0, 100)
    plt.xlabel('Time (s)')
    plt.ylabel('Smell Concentration u')
    plt.title('Smell Concentration at Professor Locations')
    plt.legend()
    plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.show()

    pass



def main():
    # question_1()
    # question_2()
    question_4()
    pass

if __name__ == "__main__":
    main()