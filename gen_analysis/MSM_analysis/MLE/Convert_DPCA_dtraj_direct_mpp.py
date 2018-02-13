
# coding: utf-8

# In[2]:

import pyemma
pyemma.__version__


# In[3]:

import os
import numpy as np
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[4]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator


# In[8]:

# get the dtraj from DPCA analysis
meta_lim='0.800'
dtraj_DPCA = np.genfromtxt('../DPCA_fromTICA-5dim/clustered_traj_mpp_'+meta_lim)
dtraj_DPCA = dtraj_DPCA.astype(int)
Nprune = np.load('../DPCA_fromTICA-5dim/Nprune.npy')
dtraj_rama = np.load('../BMSM/dtraj/dtrajs_rama_2st_allres_1D_noLHhelix.npy')[:,::Nprune]
base_dir = '/data/isilon/rudzinski/cluster_tmp/AAQAA/AAQAA_hybrid_AMBER_Go/wDB-HP_inter/NC_CA/2016_10_21/'
cc_full = np.load(base_dir+'cc_full.npy')


# In[9]:

# first cat the rama dtraj to align with the DPCA dtraj
dtraj_rama_cat = np.concatenate(dtraj_rama).astype(int)


# In[10]:

print dtraj_rama_cat.shape
print dtraj_DPCA.shape


# In[11]:

# now let's calculate some properties of each cluster


# In[12]:

# first, simply what is the fractional make-up of the cluster, in terms of the rama microstates


# In[13]:

clust_make = np.zeros(shape=(np.unique(dtraj_DPCA).shape[0],cc_full.shape[0]))
for clust_ind,clust in enumerate(np.unique(dtraj_DPCA)):
    frs = np.where( dtraj_DPCA == clust )[0]
    micros = dtraj_rama_cat[frs]
    for micro_ind, micro in enumerate(np.unique(micros)):
        clust_make[clust_ind,micro] = np.where(micros==micro)[0].shape[0] / float(len(micros))


# In[14]:

# for calculating eqm prop
state_weights_eqm_prop = np.load(base_dir+'state_weights_eqm_prop.npz')
state_weights_eqm_prop.files


# In[15]:

Nhi_weight = state_weights_eqm_prop['Nhi_weight']
Nh_weight = state_weights_eqm_prop['Nh_weight']
Nclust = clust_make.shape[0]
clust_Nh_i = np.zeros(shape=(np.unique(dtraj_DPCA).shape[0],cc_full.shape[1]))
clust_Nh = np.zeros(np.unique(dtraj_DPCA).shape[0])
for clust_ind in range(Nclust):
    clust_Nh_i[clust_ind] = np.sum(clust_make[clust_ind]*Nhi_weight.T,axis=1)
    clust_Nh[clust_ind] = np.sum(clust_make[clust_ind]*Nh_weight.T)


# In[ ]:




# In[16]:

feat = np.vstack((clust_Nh,np.arange(len(clust_Nh)))).T

clust_sorted = sorted(feat, key=lambda k: (k[0]), reverse=True)
clust_sorted = np.array(clust_sorted)
clust_map = clust_sorted[:,1].astype(int)


# In[17]:

# get the backmap from original clust number
clust_backmap_dic = {}
clust_num_OLD = np.unique(dtraj_DPCA)
for state_ind,state in enumerate(clust_map):
    clust_backmap_dic[clust_num_OLD[state]] = state_ind


# In[18]:

# reorder everything that we are planning to save
clust_make = clust_make[clust_map]


# In[20]:

# renumber the clusters and put back into traj form and save the mapping


# In[21]:

n_traj = dtraj_rama.shape[0]
traj_len = dtraj_rama.shape[1]
clust_num_OLD = np.unique(dtraj_DPCA)
dtraj_DPCA_renum = []
for traj in range(n_traj):
    dtraj_DPCA_renum.append([])
    for fr in range(traj_len):
        OG_fr = traj*traj_len + fr
        dtraj_DPCA_renum[traj].append( clust_backmap_dic[dtraj_DPCA[OG_fr]] )
    dtraj_DPCA_renum[traj] = np.squeeze(np.array(dtraj_DPCA_renum[traj]))


# In[ ]:




# In[22]:

# save all the necessities
np.save('mu_clust_mpp',clust_make)
np.save('dtraj_DPCA_mpp', dtraj_DPCA_renum)
np.save('DPCA_backmap_dic_mpp', clust_backmap_dic)


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



