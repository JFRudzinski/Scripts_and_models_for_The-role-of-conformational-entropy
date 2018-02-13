
# coding: utf-8

# AA MSM from Qhel directly via rama-plot
# -------
# with nearest neighbor correlations
# ====

# In[1]:

import pyemma
pyemma.__version__

import numpy as np

# In[2]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator
from copy import deepcopy

# Read in the opt LR params
# ------
Nres = 7
LR_opt = np.load('../../LR_opt.npz')
w_opt = np.exp(LR_opt['lnw_opt'])[0]
v_opt = np.exp(LR_opt['lnv_opt'])[0]

# convert to ZB params
sig = w_opt / (1.+v_opt)
s = (v_opt**2) / ((1.+v_opt)**4)

# and define the rates
gamma_h = 1.
gamma_c = 1.
k0 = 0.05
k1 = 1.
k2 = k1 / s
k3 = gamma_h * sig * k1
k4 = gamma_h * k2
k5 = gamma_c * sig * k2
k6 = gamma_c * k1

# get the T3 states
cc3res = np.load('../../../../DCA/cc3res.npy')
'''
[0 0 0]
[1 0 0]
[0 0 1]
[0 1 0]
[1 0 1]
[1 1 0]
[0 1 1]
[1 1 1]
'''
#for indi,si in enumerate(cc3res):
#    for indj,sj in enumerate(cc3res):
#        if ( np.sum(np.abs(si-sj)) != 1 ):
#            continue
#        print indi,indj

K3res = np.zeros(shape=(len(cc3res),len(cc3res)))
# coil nucleation: s0 [0 0 0] -> s1 [1 0 0], s2 [0 0 1], s3 [0 1 0]
K3res[0,1] = k5
K3res[0,2] = k5
K3res[0,3] = k5
K3res[1,0] = k6
K3res[2,0] = k6
K3res[3,0] = k6
# helix nucleation: s7 [1 1 1] -> s4 [1 0 1], s5 [1 1 0], s6 [0 1 1] OR s5 [1 1 0], s6 [0 1 1] -> s3 [0 1 0]
K3res[7,4] = k3
K3res[7,5] = k3
K3res[7,6] = k3
K3res[6,3] = k3
K3res[5,3] = k3
K3res[4,7] = k4
K3res[5,7] = k4
K3res[6,7] = k4
K3res[3,6] = k4
K3res[3,5] = k4
# helix elongation: s4 [1 0 1], s6 [0 1 1] -> s2 [0 0 1] OR s4 [1 0 1], s5 [1 1 0] -> s1 [1 0 0]
K3res[4,2] = k1
K3res[6,2] = k1
K3res[4,1] = k1
K3res[5,1] = k1
K3res[2,4] = k2
K3res[2,6] = k2
K3res[1,4] = k2
K3res[1,5] = k2
# now set the diagonal elements
for state in range(len(cc3res)):
    K3res[state,state] = -np.sum(K3res[state])

# Get the transition probabilities
T_zip = scipy.linalg.expm(k0*K3res)
print np.sum(T_zip,axis=1)
print np.diag(T_zip)

import matplotlib.pyplot as plt
plt.pcolor(T_zip)
plt.colorbar()
plt.show()

# Now save all the model
np.save('T3res',T_zip)

