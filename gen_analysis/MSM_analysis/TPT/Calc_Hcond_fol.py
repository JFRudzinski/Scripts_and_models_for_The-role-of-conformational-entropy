
# coding: utf-8

# In[ ]:

import pyemma
pyemma.__version__

import numpy as np

# In[ ]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})
#from matplotlib import rc
#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


# In[ ]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator
from copy import deepcopy

# Get the model
TNres_CG = np.load('../MLE/T_mle.npy')
# Get the probability that a traj passes from i to j through d
Alpha_CG = np.load('Alpha_fol.npy')

# Let's compute the conditional path entropy for each of the dominant reaction pathways
s = TNres_CG.shape[0]-1
d = 0

# get the reaction pathways from s to d
mle_Nres_CG = pyemma.msm.markov_model(TNres_CG)
tpt_mle_CG = msm.tpt(mle_Nres_CG,[s],[d])
paths_CG_fol = np.array(tpt_mle_CG.pathways())

# compute the conditional path entropy according to Kafsi et al. IEEE Inf Theory (2013)
H_path_cond_CG = []
for path in paths_CG_fol[0]:
    tmp = 1.
    H_tmp = 0.
    for node in range(len(path)-2):
        # Get Tp
        Tp = deepcopy(TNres_CG)
        for row in range(TNres_CG.shape[0]):
            if ( (row==path[node+1]) or (row==d) ):
                Tp[row] *= 0.
                Tp[row,row] = 1.
                continue
            elif ( abs(Alpha_CG[row,path[node+1]]-1.)<1e-6 ):
                continue

            for col in range(TNres_CG.shape[1]):
                Tp[row,col] *= (1.-Alpha_CG[col,path[node+1]]) / (1.-Alpha_CG[row,path[node+1]])
        if ( not np.all( np.abs(np.sum(Tp,axis=1)-1.)<1e-6) ):
            raise ValueError('Tp is not stochastic!')
        # calc H(Tp)k,k+1
        Qd = np.delete(TNres_CG,[path[node+1]],axis=0) # nb - I thought this should be Tp, but then you will always have a singular matrix...confusing notation!
        Qd = np.delete(Qd,[path[node+1]],axis=1)
        I_N = np.identity(TNres_CG.shape[0]-1)
        M_F = np.linalg.inv(I_N-Qd)
        # Get the local entropies
        S = -Tp*np.ma.log2(Tp)
        S = np.ma.MaskedArray(S,fill_value=0)
        S = np.ma.filled(S)
        S = np.sum(S,axis=1)
        m_ind = deepcopy(path[node])
        if ( m_ind == path[node+1] ): # starting and intermediate states are the same (effectively no intermediate), just move on
            H_tmp += 0.
            continue
        if ( m_ind > path[node+1] ):
            m_ind -= 1
        for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[path[node+1]]) ):
            H_tmp += M_F[m_ind,ind]*S[node2]
    # calc H(T)ld
    l = len(path)-2
    Qd = np.delete(TNres_CG,[d],axis=0) # nb - I thought this should be Tp, but then you will always have a singular matrix...confusing notation!
    Qd = np.delete(Qd,[d],axis=1)
    I_N = np.identity(TNres_CG.shape[0]-1)
    M_F = np.linalg.inv(I_N-Qd)
    # Get the local entropies and normalize
    S = -TNres_CG*np.ma.log2(TNres_CG)
    S = np.ma.MaskedArray(S,fill_value=0)
    S = np.ma.filled(S)
    S = np.sum(S,axis=1) #/ np.max(S)
    m_ind = deepcopy(path[l])
    if ( m_ind > d ):
        m_ind -= 1
    for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[d]) ):
        H_tmp += M_F[m_ind,ind]*S[node2]
    H_path_cond_CG.append(H_tmp)

np.save('H_path_cond_CG_fol.npy', H_path_cond_CG)



# repeat the calculation for single intermediate states; i.e., H_{sd|u} for each u
# this is the entropy of reactive trajectories given that they pass through u
H_path_cond_state_CG = np.zeros(TNres_CG.shape[0])
for state in range(1,TNres_CG.shape[0]-1):
    path = np.array([s,state,d])
    tmp = 1.
    H_tmp = 0.
    for node in range(len(path)-2):
        # Get Tp
        Tp = deepcopy(TNres_CG)
        for row in range(TNres_CG.shape[0]):
            if ( (row==path[node+1]) or (row==d) ):
                Tp[row] *= 0.
                Tp[row,row] = 1.
                continue
            elif ( abs(Alpha_CG[row,path[node+1]]-1.)<1e-6 ):
                continue

            for col in range(TNres_CG.shape[1]):
                Tp[row,col] *= (1.-Alpha_CG[col,path[node+1]]) / (1.-Alpha_CG[row,path[node+1]])
        if ( not np.all( np.abs(np.sum(Tp,axis=1)-1.)<1e-6) ):
            raise ValueError('Tp is not stochastic!')
        # calc H(Tp)k,k+1
        Qd = np.delete(TNres_CG,[path[node+1]],axis=0) # nb - I thought this should be Tp, but then you will always have a singular matrix...confusing notation!
        Qd = np.delete(Qd,[path[node+1]],axis=1)
        I_N = np.identity(TNres_CG.shape[0]-1)
        M_F = np.linalg.inv(I_N-Qd)
        # Get the local entropies
        S = -Tp*np.ma.log2(Tp)
        S = np.ma.MaskedArray(S,fill_value=0)
        S = np.ma.filled(S)
        S = np.sum(S,axis=1)
        m_ind = deepcopy(path[node])
        if ( m_ind == path[node+1] ): # starting and intermediate states are the same (effectively no intermediate), just move on
            H_tmp += 0.
            continue
        if ( m_ind > path[node+1] ):
            m_ind -= 1
        for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[path[node+1]]) ):
            H_tmp += M_F[m_ind,ind]*S[node2]
    # calc H(T)ld
    l = len(path)-2
    Qd = np.delete(TNres_CG,[d],axis=0) # nb - I thought this should be Tp, but then you will always have a singular matrix...confusing notation!
    Qd = np.delete(Qd,[d],axis=1)
    I_N = np.identity(TNres_CG.shape[0]-1)
    M_F = np.linalg.inv(I_N-Qd)
    # Get the local entropies and normalize
    S = -TNres_CG*np.ma.log2(TNres_CG)
    S = np.ma.MaskedArray(S,fill_value=0)
    S = np.ma.filled(S)
    S = np.sum(S,axis=1) #/ np.max(S)
    m_ind = deepcopy(path[l])
    if ( m_ind > d ):
        m_ind -= 1
    for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[d]) ):
        H_tmp += M_F[m_ind,ind]*S[node2]
    H_path_cond_state_CG[state] = 1.*H_tmp

np.save('H_path_cond_state_CG_fol.npy', H_path_cond_state_CG)


