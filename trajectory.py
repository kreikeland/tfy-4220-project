import numpy as np
from magnetic import B, M
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import argparse
from scipy.integrate import ode, odeint, solve_ivp


R_e = 6.37e6 # radius of Earth, m
e = 1.602e-19 # elementary charge, C
m_p = 1.67e-27 # proton mass, kg
c = 3e8 # speed of light, units?

def compute_trajectory_rk4(x0, v0, t0, tf, dt, q, m):
    """
    3D Runge-Kutta 4th order solver for particle trajectory.
    Inspired by https://arxiv.org/pdf/1112.3487.pdf
    
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
    vsq = np.linalg.norm(v0)**2
    gamma = 1/np.sqrt(1-vsq/c**2) # gamma approximately constant for entire trajectory
    
    def rhs(t, Y):
       x,y,z,vx,vy,vz = Y
       Bx, By, Bz = B(x,y,z, M)
       fac = q/(m*gamma)
       return [ vx, vy, vz,
                    fac*(vy*Bz - vz*By),
                    fac*(vz*Bx - vx*Bz),
                    fac*(vx*By - vy*Bx) ]
    
    res = solve_ivp(rhs, [t0, tf], np.array([x0[0], x0[1], x0[2], v0[0], v0[1], v0[2]]))

    return res.t.T, res.y.T # transpose to get Nx6 array


def plot_traj_3d(r, K):
    '''
    Plot the 3d trajectory given the position vectors r

    input:
        r: 3xNsteps array

    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Trajectory for {K*1e-6}MeV proton with $x_0 = {r[0]}R_E$')
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
    K = 1e7 # kinetic energy, eV
    Kj = K*e   # convert to Joule
    v_mag = c*np.sqrt(1-(m_p*c**2)**2/(m_p*c**2 + Kj)**2)
	
    # pitch_angle = 30.0 # initial angle between velocity and mag.field (degrees)
    # vz0 = v_mag*np.cos(pitch_angle*np.pi/180)
    # vy0 = v_mag*np.sin(pitch_angle*np.pi/180)
    # vx0 = 0
    # v0 = np.array([vx0,vy0,vz0])
    v0 = np.array([1,0,0])*v_mag
    t0 = 0
    tf = 120
    dt = 0.01
    x0 = np.array([-2.0,0,0])*R_e 
    q = e
    m = m_p

    t, output = compute_trajectory_rk4(x0, v0, t0, tf, dt, q, m)
    r, v = output[:,0:3], output[:,3:6]
    plot_traj_3d(r/R_e, K)
    

