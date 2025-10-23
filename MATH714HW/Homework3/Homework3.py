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


def main():
    question_1()

if __name__ == "__main__":
    main()