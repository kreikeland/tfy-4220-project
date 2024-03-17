import numpy as np
import os

nx, ny, nz = 11, 3, 3
x, y, z = np.linspace(-2, -12, nx).round(1), np.linspace(-1, 1, ny).round(1), np.linspace(-1, 1, nz).round(1)

with open('initial_conditions/initial_conditions.txt', 'w') as file:
    file.write('x0 v0\n')
    for i in range(ny):
        for j in range(nz):
            for k in range(nx):
                x0 = f'{x[k]},{y[i]},{z[j]}'
                v0 = f'1,0,0'
                file.write(f'{x0} {v0}\n')
    file.close()