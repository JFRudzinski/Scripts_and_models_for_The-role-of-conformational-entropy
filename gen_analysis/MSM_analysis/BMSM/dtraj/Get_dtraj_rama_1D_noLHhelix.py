
# coding: utf-8

# CG MSM from Qhel directly via rama-plot
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

import os
import numpy as np
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


# Read in the dtrajs
# ------

# In[4]:

dtraj_rama = np.load('traj_rama.npy')


# In[5]:

Nrama = dtraj_rama.shape[0]
Ntraj = dtraj_rama.shape[1]
Nfr = dtraj_rama.shape[2]
Ndih = dtraj_rama.shape[3]


# In[7]:

dtraj_rama.shape


# In[4]:




# In[8]:

for rama in range(Nrama):
    for i in range( Ntraj ):
        dtraj_rama[rama][i][np.where(dtraj_rama[rama][i][:,1] < -125)[0],1] += 360  


# In[9]:

dtraj_phi = []
for rama in range(Nrama):
    dtraj_phi.append([])
    for i in range( len(dtraj_rama[rama]) ):
        dtraj_phi[rama].append(dtraj_rama[rama][i][:,1])


# **simple clustering along psi only for discretization**

# In[10]:

n_clusters = 2     # number of k-means clusters


# In[11]:

clustering_rama = []
for rama in range(Nrama):
    clustering_rama.append(coor.cluster_kmeans(dtraj_phi[rama],k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True))

dtrajs_rama = []
for rama in range(Nrama):
    dtrajs_rama.append(clustering_rama[rama].dtrajs)


# In[14]:

for rama in range(Nrama):
    for traj in range( len(dtraj_rama[rama]) ):
        if ( dtrajs_rama[rama][traj][np.where(dtraj_phi[rama][traj] < 0)[0][0]] != 0 ):
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == 0)[0]] -= 1
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == 1)[0]] -= 1
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == -1)[0]] += 2

dtraj_psi = []
for rama in range(Nrama):
    dtraj_psi.append([])
    for i in range( len(dtraj_rama[rama]) ):
        dtraj_psi[rama].append(dtraj_rama[rama][i][:,0])


# remove the LH-helix states from the helix
for rama in range(Nrama):
    for traj in range( Ntraj ):
        LH_frs = np.where( dtraj_psi[rama][traj] > 0 )[0]
        dtrajs_rama[rama][traj][LH_frs] = 1

dtrajs = []
for traj in range( len(dtrajs_rama[0]) ):
    tmp = dtrajs_rama[0][traj]
    for rama in range(1,len(dtrajs_rama)):
        tmp = np.vstack( (tmp,dtrajs_rama[rama][traj]) )
    tmp = tmp.T
    dtrajs.append( tmp )
    dtrajs[traj].astype('int64')


# In[23]:

dtrajs[0].shape


# In[24]:

dtrajs_sum = []
for i in range( len(dtrajs) ):
    dtrajs_sum.append( np.sum(dtrajs[i],axis=1) )
    dtrajs_sum[i].astype('int64')

# we need a single dimensional identifier of the microstate, can we cluster to automize?
n_clusters = 16
clustering = coor.cluster_regspace(dtrajs_sum,max_centers=n_clusters,dmin=0.5)
#clustering = coor.cluster_kmeans(dtrajs,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)


# In[40]:

dtrajs_1D = clustering.dtrajs


# In[41]:

cc = clustering.clustercenters[:]
cc


# In[42]:

sorted_list = sorted(cc)
cc_sorted = np.array(sorted_list)[:]
cc_sorted = cc_sorted.tolist()
for i in range(len(cc)):
    cc_sorted[i] = map(int,cc_sorted[i])
cc_sorted


# In[43]:

for traj in range(len(dtrajs_1D)):
    dtrajs_1D[traj] = cc[dtrajs_1D[traj]].astype(int)


# In[25]:

# we need a single dimensional identifier of the microstate, can we cluster to automize?
# n_clusters = 2**15
# clustering = coor.cluster_regspace(dtrajs,max_centers=n_clusters,dmin=0.5)
#clustering = coor.cluster_kmeans(dtrajs,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)

# this is unnecessary and expensive for so many states, we already generated the states, just read them in
cc = np.load('/data/isilon/rudzinski/cluster_tmp/AAQAA/AAQAA_hybrid_AMBER_Go/wDB-HP_inter/NC_CA/2016_10_21/cc_full.npy')


# In[26]:

#dtrajs_full_1D = clustering.dtrajs


# In[27]:

# cc = clustering.clustercenters[:]
cc.shape


# In[28]:

# convert the full traj to 1D
dtrajs_full_1D = []
for traj in range(len(dtrajs)):
    dtrajs_full_1D.append([])
    for fr in range(dtrajs[traj].shape[0]):
        state = np.where( np.all(dtrajs[traj][fr]==np.array(cc),axis=1) == True )[0][0]
        dtrajs_full_1D[traj].append(state)
    dtrajs_full_1D[traj] = np.array(dtrajs_full_1D[traj])
    print 'done with traj '+str(traj)


np.save('dtrajs_rama_2st_allres_1D_noLHhelix', dtrajs_full_1D)
np.save('dtrajs_Qhel_allres_1D_noLHhelix', dtrajs_1D)


