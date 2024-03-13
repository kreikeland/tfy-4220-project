import numpy as np
from magnetic import B, M
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import argparse

R_e = 6.37e6 # radius of Earth, 10⁶m
e = 1.602e-19 # elementary charge, C
m_p = 1.67e-27 # proton mass, kg
c = 10e8/R_e # speed of light, units?

def compute_trajectory_rk4(x0, v0, t0, tf, dt, q, m):
    """
    3D Runge-Kutta 4th order solver for particle trajectory.
    
    Parameters:
        x0: array-like
            The initial position vector.
        v0: array-like
            The initial velocity vector
        t0: float
            The initial time.
        tf: float
            The final time.
        dt: float
            The time step size.
        q: float
            Particle charge
        m: float
            Particle mass
    
    Returns:
        t: array-like
            The array of time values.
        x: array-like
            The array of position vectors at each time step.
        v: array-like
            The array of velocity vectors at each time step.
    
    """    
    t = np.arange(t0, tf + dt, dt)
    n = len(t)
    x = np.zeros((n, len(x0)))
    x[0] = x0
    v = np.zeros((n, len(v0)))
    v[0] = v0

    for i in range(n - 1):
            k1_v = dt * rhs_v(x[i], v[i], q, m)
            k2_v = dt * rhs_v(x[i], v[i] + 0.5 * k1_v, q, m)
            k3_v = dt * rhs_v(x[i], v[i] + 0.5 * k2_v, q, m)
            k4_v = dt * rhs_v(x[i], v[i] + k3_v, q, m)
            v[i + 1] = v[i] + (1/6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

            k1_r = dt * rhs_r(v[i])
            k2_r = dt * rhs_r(v[i] + 0.5 * k1_r)
            k3_r = dt * rhs_r(v[i] + 0.5 * k2_r)
            k4_r = dt * rhs_r(v[i] + k3_r)
            x[i + 1] = x[i] + (1/6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)

    return t, x, v

def work(q,r,v):
    '''
    Calculate the work done by the magnetic force for a given trajectory
    r, v
    Test to gauge the accuracy of numerical method
    '''
    Bvec = np.array(B(r[:,0], r[:,1], r[:,2], M)).T
    force = q*(np.cross(v, Bvec))
    return (force*v).sum(1)

def rhs_v(x, v, q, m):
    '''
    RHS of EOM for v
    '''
    gamma = 1 / (np.sqrt(1-np.linalg.norm(v)**2/c**2))
    Bvec = np.array(B(x[0], x[1], x[2], M))
    return q / (gamma*m)* np.cross(v, Bvec) 

def rhs_r(v):
    '''
    RHS of EOM for r
    '''
    return v

def plot_traj_3d(r):
    '''
    Plot the 3d trajectory given the position vectors r

    input:
        r: 3xNsteps array

    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(r[:,0], r[:,1], r[:,2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    # draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color="r")
    plt.show()


if __name__ == '__main__':
    # assume unit mass and charge
    # beware of units! magnetic field calculated in length units of R_e
    
    v0 = np.array([1/np.sqrt(2),1/np.sqrt(2),0])*0.05 # units? R_e s⁻1
    t0 = 0
    tf = 10000
    dt = 0.1
    x0 = np.array([2,0,0]) # units? R_e?
    q = 1e3
    m = 1

    t, r, v = compute_trajectory_rk4(x0, v0, t0, tf, dt, q, m)
    w = work(q,r,v)
    print(np.sum(w))
    

