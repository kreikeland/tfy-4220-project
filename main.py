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
            x0 = np.array(x0.split(','), dtype=float)
            v0 = np.array(v0.split(','), dtype=float)
            init_conditions.append((x0, v0))
    return init_conditions

def run_sim(initial_conditions, v_mag, tf):
    '''
    Compute trajectories of particles with initial velocity v_mag 
    input:
        initial_conditions:
            N x 2 x 3 array of initial conditions. First entry x0, second v0
        v_mag: 
            float. initial velocity magnitude for all runs.
        tf:
            float. final time of simulation. 
    '''
    dt = 0.1
    t0 = 0
    q = 1000
    m = 1
    with alive_bar(len(initial_conditions)) as bar:
        for i, x0v0 in enumerate(initial_conditions):
            x0, v0 = x0v0[0], x0v0[1]
            # run computation
            t, x, v = compute_trajectory_rk4(x0, v0 * v_mag, t0, tf, dt, q, m)
            # save trajectories
            savepath = os.path.join(os.getcwd(), 'trajectories/')
            np.savetxt(savepath+f'trajectory_r{i}.txt', x, delimiter=',')
            np.savetxt(savepath+f'trajectory_v{i}.txt', v, delimiter=',')
            bar()
      
if __name__ == "__main__":
    init_conditions = read_initial_conditions('initial_conditions.txt')
    run_sim(init_conditions, 0.05, 5000)
    # x = np.loadtxt('trajectories/trajectory_r0.txt', delimiter=',')
    # plot_traj_3d(x)