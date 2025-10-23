import numpy as np 
import matplotlib.pyplot as plt

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
    # # We need to plot the stability region for the TR-BDF2 method where z(\theta) &= \frac{7e^{i\theta} + 5 \pm \sqrt{(7e^{i\theta} + 5)^2 - 48(e^{i\theta} - 1)}}{2e^{i\theta}}.
    # The plot the following code plots is incorrect.
    theta = np.linspace(0, 2*np.pi, 400)
    z_plus = (7*np.exp(1j*theta) + 5 + np.sqrt((7*np.exp(1j*theta) + 5)**2 - 48*(np.exp(1j*theta) - 1))) / (2*np.exp(1j*theta))
    z_minus = (7*np.exp(1j*theta) + 5 - np.sqrt((7*np.exp(1j*theta) + 5)**2 - 48*(np.exp(1j*theta) - 1))) / (2*np.exp(1j*theta))

    plt.figure(figsize=(6,6))
    plt.plot(z_plus.real, z_plus.imag, 'b-', label='Branch +')
    plt.plot(z_minus.real, z_minus.imag, 'r-', label='Branch -')
    plt.fill_between(np.concatenate([z_plus.real, z_minus.real[::-1]]),
                     np.concatenate([z_plus.imag, z_minus.imag[::-1]]),
                     color='lightgray', alpha=0.5)
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(0, color='black', lw=0.5)
    plt.title('Region of Absolute Stability for the TR-BDF2 Method')
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    plt.show()

    # Code below can be uncommented to plot |R(z)| <= 1 instead of using the boundary locus method

    # Plot |R(z)| <= 1
    def R(z):
        return (5*z + 12) / ((z - 4)*(z - 3))
    
    # Plotting on a -20 to 20 real and -15 to 15 imaginary grid

    x = np.linspace(-20, 20, 400)
    y = np.linspace(-20, 20, 400)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j*Y
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



def main():
    # question_1()
    question_2()
    pass

if __name__ == "__main__":
    main()