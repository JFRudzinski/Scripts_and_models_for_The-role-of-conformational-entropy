import numpy as np

Ntraj = 40
traj_len = np.genfromtxt('traj_len.dat')
Nprune = np.load('Nprune.npy')

pop_min = (0.0001*Ntraj*traj_len / Nprune).astype(int)

print pop_min

np.savetxt('pop_min.dat',[pop_min],fmt='%d')
