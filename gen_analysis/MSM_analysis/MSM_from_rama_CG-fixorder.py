
# coding: utf-8

# MSM from "primitive clustering" along rama-plot
# ====

# In[2]:

import pyemma
pyemma.__version__


# In[3]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})

import numpy as np

# In[4]:

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

# In[10]:

indir = '/data/isilon/rudzinski/cluster_tmp/AAQAA/AAQAA_hybrid_AMBER_Go/wDB-HP_inter/NC_CA/2016_10_21/epsNC-11/epsNC-11_epsdb-0.3epsNC_epshp-0.35epsNC/T-300/'

traj_dir_base = 'run_from_Qhel-'
Qhel_val = ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9']

Nrama = 15
GLN_rama = [4,9,14]
dtraj_rama = []
for rama in range(Nrama):
    dtraj_rama.append([])
    for i in range (len(Qhel_val)):
        for j in range(4):
            traj_dir = indir+traj_dir_base+Qhel_val[i]+'/run'+str(j+1)+'/'
            if (rama+2 in GLN_rama):
                dtraj_rama[rama].append( np.genfromtxt(traj_dir+'rama_GLN'+str(rama+2)+'.xvg')[:,0:2] )
            else:
                dtraj_rama[rama].append( np.genfromtxt(traj_dir+'rama_ALA'+str(rama+2)+'.xvg')[:,0:2] )


# In[11]:

for rama in range(Nrama):
    for i in range( len(dtraj_rama[0]) ):
        dtraj_rama[rama][i][np.where(dtraj_rama[rama][i][:,1] < -125)[0],1] += 360  


# In[12]:

dtraj_phi = []
for rama in range(Nrama):
    dtraj_phi.append([])
    for i in range( len(dtraj_rama[rama]) ):
        dtraj_phi[rama].append(dtraj_rama[rama][i][:,1])


# **simple clustering along psi only for discretization**

# In[13]:

n_clusters = 2     # number of k-means clusters


# In[14]:

clustering_rama = []
for rama in range(Nrama):
    clustering_rama.append(coor.cluster_kmeans(dtraj_phi[rama],k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True))


# In[15]:

cc_rama = []
for rama in range(Nrama):
    cc_rama.append(clustering_rama[rama].clustercenters[:,0])
    print cc_rama[rama]


# In[16]:

dtrajs_rama = []
for rama in range(Nrama):
    dtrajs_rama.append(clustering_rama[rama].dtrajs)


# In[17]:

for rama in range(Nrama):
    for traj in range( len(dtraj_rama[rama]) ):
        if ( dtrajs_rama[rama][traj][np.where(dtraj_phi[rama][traj] < 0)[0][0]] != 0 ):
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == 0)[0]] -= 1
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == 1)[0]] -= 1
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == -1)[0]] += 2


# In[18]:

dtrajs = []
for i in range( len(dtraj_rama[0]) ):
    tmp = np.stack( (dtrajs_rama[0][i], dtrajs_rama[1][i]) )
    for rama in range(2,Nrama):
        tmp = np.vstack( (tmp, dtrajs_rama[rama][i]) )
    dtrajs.append( tmp.T )
    dtrajs[i].astype('int64')


# In[24]:

Nh = 1.0 - np.mean(dtrajs)
print 'Nh ='+str(Nh)


# In[29]:

from copy import deepcopy
dtrajs_LR = deepcopy(dtrajs)


# In[31]:

for traj in range(len(dtrajs)):
    for fr in range(dtrajs[traj].shape[0]):
        for rama in range(dtrajs[traj].shape[1]):
            if ( dtrajs[traj][fr][rama] == 1 ):
                dtrajs_LR[traj][fr][rama] = 0
            else:
                dtrajs_LR[traj][fr][rama] = 1
                if ( (rama!=0) and (rama!=dtrajs[traj].shape[1]-1) ):
                    if ( (dtrajs[traj][fr][rama-1]==0) and (dtrajs[traj][fr][rama+1]==0) ):
                        dtrajs_LR[traj][fr][rama] = 2


# In[40]:

N_h_LR = 0
for traj in range(len(dtrajs)):
    N_h_LR += len(np.where(dtrajs_LR[traj]==2)[0])
N_h_LR /= 1.*len(dtrajs)*dtrajs[0].shape[0]*(dtrajs[0].shape[1]-2)


# In[41]:

print '\n'
print 'N_h_LR = '+str(N_h_LR)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



