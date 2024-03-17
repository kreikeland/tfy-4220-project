import numpy as np
import argparse
from trajectory import plot_traj_3d, R_e
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from initial_conditions.init_cond import nx, ny, nz
from alive_progress import alive_bar

# generates 2d plots of trajectories for 100eV, 10keV, 1MeV
if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-sim', default=3)
    # parser.add_argument('-id', default=0)
    # args = parser.parse_args()
    # sim = args.sim
    # id = args.id
    selection = [36]
    with alive_bar(nx*ny*nz) as bar:
        # for i in range(nx*ny*nz):  # loop over all
        for i in selection:         # loop over the selected, bar acts weird    
            traj = np.loadtxt(f'trajectories/trajectory_sim{1}_K-1000000.0eV_{i}.txt', delimiter=',')
            t, x, v = traj[:,0], traj[:,1:4], traj[:,4:7]
            
            traj2 = np.loadtxt(f'trajectories/trajectory_sim{1}_K-10000.0eV_{i}.txt', delimiter=',')
            t2, x2, v2 = traj2[:,0], traj2[:,1:4], traj2[:,4:7]

            traj3 = np.loadtxt(f'trajectories/trajectory_sim{1}_K-100.0eV_{i}.txt', delimiter=',')
            t3, x3, v3 = traj3[:,0], traj3[:,1:4], traj3[:,4:7]

            fig, ax = plt.subplots()
            ax.plot(x[:,0]/R_e, x[:,1]/R_e, label='1MeV')
            ax.plot(x2[:,0]/R_e, x2[:,1]/R_e, label='10keV')
            ax.plot(x3[:,0]/R_e, x3[:,1]/R_e, label='100eV')
            ax.add_patch(Circle((0, 0), 1, color='r', zorder=0))
            ax.set_xlabel('x [R_e]')
            ax.set_ylabel('y [R_e]')
            ax.legend()
            print("Initial position: ",x[0]/R_e, "R_e", "id: ", i)
            plt.savefig(f'figures/trajectories-xy-{i}.eps', bbox_inches='tight')
            plt.close(fig)
            bar()