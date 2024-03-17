import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
import matplotlib.cm as cm
import argparse

# constants
mu0 = 4*np.pi*1e-7 # permeability of free space, T*m/A
B0 = 3.07e-5 # magnetic field at Earth's surface, T
R_e = 6.37e6 # radius of Earth, m
alpha = (23+11) * np.pi / 180 # angle of inclination of Earth's magnetic field, 11deg + 23deg
# alpha=0
M = R_e**3 * B0 * np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                        [0,1,0],
                        [np.sin(alpha), 0, np.cos(alpha)]]) @ np.array([0,0,1]) # Earth's magnetic dipole moment, T x m**3, rotated

def B(x,y,z, m):
    """
    Return the Earths magnetic field components at a point xyz
    input:
        x, y, z: float. position coordinates
        m: array-like. magnetic dipole moment
    return:
        Bx, By, Bz: float. magnetic field components    
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    Bx = (3*x**2 - r**2) * m[0] + 3*x*y * m[1] + 3*x*z * m[2]
    By = 3*x*y * m[0] + (3*y**2 - r**2) * m[1] + 3*y*z * m[2]
    Bz = 3*x*z * m[0] + 3*y*z * m[1] + (3*z**2 - r**2) * m[2]
    return Bx/r**5, By/r**5, Bz/r**5


if __name__ == '__main__':
    # Plot xz-plane
    x, z = np.linspace(-6*R_e, 6*R_e, 100), np.linspace(-6*R_e, 6*R_e, 100)
    X, Z = np.meshgrid(x, z)
    Y = 0
    Bx, By, Bz = B(X, Y, Z, M)
    B_ = np.sqrt(Bx**2 + By**2 +  Bz**2)

    fig, ax = plt.subplots()
    ax.streamplot(X/R_e, Z/R_e, Bx, Bz, color='k',linewidth=1,
                density=2.5, arrowstyle='->', arrowsize=1)
    ax.add_patch(Circle((0, 0), 1, color='k', zorder=100))
    # p = ax.contourf(X, Z, By, cmap=cm.jet, levels=0) # contour for By, pos/neg
    # c = fig.colorbar(p)
    # c.set_label('$B_{y}$')
    ax.set_aspect('equal')
    ax.set_xlabel('$x$ [$R_{E}$]')
    ax.set_ylabel('$z$ [$R_{E}$]')
    ax.text(-8,6,'(a)',fontsize=12)

    #ax.set_title(f'Magnetic field lines for y={Y}$R_E$')
    fig.savefig('figures/magnetic-field-lines-x0z.eps', bbox_inches='tight')
    # plt.show()

    # Plot xy-plane
    x, y = np.linspace(-6*R_e, 6*R_e, 100), np.linspace(-6*R_e, 6*R_e, 100)
    X, Y = np.meshgrid(x, y)
    Z = -1*R_e
    Bx, By, Bz = B(X, Y, Z, M)

    fig, ax = plt.subplots()
    ax.streamplot(X/R_e, Y/R_e, Bx, By, color='k',linewidth=1,
                density=2.5, arrowstyle='->', arrowsize=1)
    # ax.add_patch(Circle((0, 0), 1, color='b', zorder=100))
    # p = ax.contourf(X, Y, Bz, cmap=cm.jet, levels=0) # contour map for Bz, pos/neg
    # c = fig.colorbar(p)
    # c.set_label('$B_{z}$')
    ax.set_aspect('equal')
    ax.set_xlabel('$x$ [$R_{E}$]')
    ax.set_ylabel('$y$ [$R_{E}$]')
    ax.text(-8,6,'(b)', fontsize=12)
    #ax.set_title(f'Magnetic field lines xy-plane for z={np.round(Z/R_e)}$R_E$')
    fig.savefig('figures/magnetic-field-lines-xy-1.eps', bbox_inches='tight')
    # plt.show()