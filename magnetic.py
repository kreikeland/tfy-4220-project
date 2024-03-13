import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
import matplotlib.cm as cm
import argparse

# constants
alpha = (23+11) * np.pi / 180 # angle of inclination of Earth's magnetic field, 11deg + 23deg
M = 30.4e-6 * np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                        [0,1,0],
                        [np.sin(alpha), 0, np.cos(alpha)]]) @ np.array([0,0,1]) # Earth's magnetic dipole moment, µT * R_e**3, rotated

def B(x,y,z, m):
    """
    Return the Earths magnetic field components at a point xyz
    from http://www.ss.ncu.edu.tw/~yhyang/104-1/Ref-KivelsonRussell-p165-167.pdf
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    Bx = (3*x**2 - r**2) * m[0] + 3*x*y * m[1] + 3*x*z * m[2]
    By = 3*x*y * m[0] + (3*y**2 - r**2) * m[1] + 3*y*z * m[2]
    Bz = 3*x*z * m[0] + 3*y*z * m[1] + (3*z**2 - r**2) * m[2]
    return Bx, By, Bz


if __name__ == '__main__':
    # Plot xz-plane
    x, z = np.linspace(-6, 6, 100), np.linspace(-6, 6, 100)
    X, Z = np.meshgrid(x, z)
    Y = 0
    Bx, By, Bz = B(X, Y, Z, M)


    fig, ax = plt.subplots()
    ax.streamplot(X, Z, Bx, Bz, color='k',linewidth=1,
                density=2, arrowstyle='->', arrowsize=1.5)
    ax.add_patch(Circle((0, 0), 1, color='k', zorder=100))
    # p = ax.contourf(X, Z, By, cmap=cm.jet, levels=0) # contour for By, pos/neg
    # c = fig.colorbar(p)
    # c.set_label('$B_{y}$')
    ax.set_aspect('equal')
    ax.set_xlabel('$x$ [$R_{E}$]')
    ax.set_ylabel('$z$ [$R_{E}$]')
    ax.set_title(f'Magnetic field lines for y={Y}$R_E$')
    plt.show()

    # Plot xy-plane
    x, y = np.linspace(-6, 6, 100), np.linspace(-6, 6, 100)
    X, Y = np.meshgrid(x, y)
    Z = 0
    Bx, By, Bz = B(X, Y, Z, M)

    fig, ax = plt.subplots()
    ax.streamplot(X, Y, Bx, By, color='k',linewidth=1,
                density=2, arrowstyle='->', arrowsize=1.5)
    ax.add_patch(Circle((0, 0), 1, color='k', zorder=100))
    # p = ax.contourf(X, Y, Bz, cmap=cm.jet, levels=0) # contour map for Bz, pos/neg
    # c = fig.colorbar(p)
    # c.set_label('$B_{z}$')
    ax.set_aspect('equal')
    ax.set_xlabel('$x$ [$R_{E}$]')
    ax.set_ylabel('$y$ [$R_{E}$]')
    ax.set_title(f'Magnetic field lines xy-plane for z={Z}$R_E$')
    plt.show()