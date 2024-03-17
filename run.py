from trajectory import *
import numpy as np
from alive_progress import alive_bar
import sys, os


def read_initial_conditions(filename):
    init_conditions = []
    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        for line in file:
            x0, v0 = line.strip().split()
            x0 = np.array(x0.split(','), dtype=float) * R_e
            v0 = np.array(v0.split(','), dtype=float)
            init_conditions.append((x0, v0))
    return init_conditions

def run_sim(initial_conditions, K, tf, sim):
    '''
    Compute trajectories of particles with initial velocity v_mag 
    input:
        initial_conditions:
            N x 2 x 3 array of initial conditions. First entry x0, second v0
        K: 
            float. Kinetic energy.
        tf:
            float. final time of simulation. 
        sim:
            string. simulation number corresponding to initial conditions number
    '''
    dt = 0.01
    t0 = 0
    q = e
    m = m_p
    Kj = K*e   # convert to Joule
    v_mag = c*np.sqrt(1-(m*c**2)**2/(m*c**2 + Kj)**2)
    with alive_bar(len(initial_conditions)) as bar: # progress bar
        for i, x0v0 in enumerate(initial_conditions):
    
            x0, v0 = x0v0[0], x0v0[1]
            
            # run computation
            t, output = compute_trajectory_rk4(x0, v0 * v_mag, t0, tf, dt, q, m)
            x, v = output[:,0:3], output[:,3:6]
            # save trajectories
            savepath = os.path.join(os.getcwd(), 'trajectories/')
            traj = np.concatenate((t[:,np.newaxis], x, v), axis=1)  
            np.savetxt(savepath+f'trajectory_sim{sim}_K-{K}eV_{i}.txt', traj, delimiter=',') 
            # np.savetxt(savepath+f'trajectory_r_sim{sim}_K-{K}eV_{i}.txt', x, delimiter=',')
            bar()

# Run this script to compute trajectories of particles with initial conditions
# given in initial_conditions.txt and energy K.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-sim', default=0)
    args = parser.parse_args()

    K = 1e2 # kinetic energy, eV
    tf = 120 # final time of simulation, s

    # pass argument -sim !=0 to run sim
    # args.sim is passed to run_sim to save trajectories with sim number
    if args.sim != 0:
        init_conditions = read_initial_conditions(f'initial_conditions/initial_conditions.txt') # read initial conditions
        run_sim(init_conditions, K, tf, args.sim) # run sim for 10MeV protons for 120s
