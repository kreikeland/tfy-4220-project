import numpy as np
import matplotlib.pyplot as plt
from trajectory import m_p, c, R_e, e
from initial_conditions.init_cond import nx, ny, nz

def K_traj(v, m):
    '''
    Compute the kinetic energy of a particle with mass m and velocity v
    return:
        K, eV
    '''
    gamma = 1/np.sqrt(1-np.linalg.norm(v, axis=1)**2 / c**2)
    return m * c**2 * (gamma - 1) / e

if __name__ == '__main__':

    ids = [44,46,48,50,52,54]
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(9,9))
    
    for i in ids:    
        traj  = np.loadtxt(f'trajectories/trajectory_sim1_K-1000000.0eV_{i}.txt', delimiter=',')
        t, x, v = traj[:,0], traj[:,1:4], traj[:,4:7]
        K = K_traj(v, m_p)
        ax1.plot(t, np.log(K*1e-6), label=f'$x_0=${x[0]/R_e}$R_E$')
        
        traj2  = np.loadtxt(f'trajectories/trajectory_sim1_K-10000.0eV_{i}.txt', delimiter=',')
        t2, x2, v2 = traj2[:,0], traj2[:,1:4], traj2[:,4:7]
        K2 = K_traj(v2, m_p)
        ax2.plot(t2, np.log(K2*1e-4), label=f'$x_0=${x2[0]/R_e}$R_E$')
        
        traj3  = np.loadtxt(f'trajectories/trajectory_sim1_K-100.0eV_{i}.txt', delimiter=',')
        t3, x3, v3 = traj3[:,0], traj3[:,1:4], traj3[:,4:7]
        K3 = K_traj(v3, m_p)
        ax3.plot(t3, np.log(K3*1e-2), label=f'$x_0=${x3[0]/R_e}$R_E$')
    
    ax2.sharex(ax1)
    ax3.sharex(ax1)
    # ax1.sharey    ax3)
    # ax2.sharey    ax3)
    ax3.set_xlabel('t [s]')
    ax1.set_ylabel('$log$ K/1MeV')
    ax2.set_ylabel('$log$ K/10keV')
    ax3.set_ylabel('$log$ K/100eV')
    ax1.text(-0.08, 1, '(a)', fontsize=12,horizontalalignment='center', 
             verticalalignment='center', transform=ax1.transAxes)
    ax2.text(-0.08, 1, '(b)', fontsize=12,horizontalalignment='center', 
             verticalalignment='center', transform=ax2.transAxes)
    ax3.text(-0.08, 1, '(c)', fontsize=12,horizontalalignment='center', 
             verticalalignment='center', transform=ax3.transAxes)
    
    ax1.legend()
    
    # fig.tight_layout()
    plt.savefig('figures/energy.eps', bbox_inches='tight')
    # plt.show()