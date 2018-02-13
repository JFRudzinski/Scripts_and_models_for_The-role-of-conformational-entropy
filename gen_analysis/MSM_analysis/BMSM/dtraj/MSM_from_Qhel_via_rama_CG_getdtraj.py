
# coding: utf-8

# CG MSM from Qhel directly via rama-plot
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


# Read in the dtrajs
# ------

# In[4]:
# get the indir from the current dir
cdir = os.getcwd()
flag = False
wdir = 'MSM_analysis'
dir_len = len(wdir)
ctr = 0
while (not flag):
    if ( cdir[ctr:dir_len+ctr] == wdir ):
        flag = True
        dir_ind = ctr
    ctr += 1

indir = cdir[:ctr+dir_len-(len(wdir)+1)]

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


# In[5]:

for rama in range(Nrama):
    for i in range( len(dtraj_rama[0]) ):
        dtraj_rama[rama][i][np.where(dtraj_rama[rama][i][:,1] < -125)[0],1] += 360  


# In[6]:

dtraj_phi = []
for rama in range(Nrama):
    dtraj_phi.append([])
    for i in range( len(dtraj_rama[rama]) ):
        dtraj_phi[rama].append(dtraj_rama[rama][i][:,1])


# **simple clustering along psi only for discretization**

# In[7]:

n_clusters = 2     # number of k-means clusters


# In[8]:

clustering_rama = []
for rama in range(Nrama):
    clustering_rama.append(coor.cluster_kmeans(dtraj_phi[rama],k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True))


# In[9]:

cc_rama = []
for rama in range(Nrama):
    cc_rama.append(clustering_rama[rama].clustercenters[:,0])
    #plt.plot(cc_rama[rama],np.zeros(len(cc_rama[rama])),marker='x')
    #print cc_rama[rama]


# In[11]:

dtrajs_rama = []
for rama in range(Nrama):
    dtrajs_rama.append(clustering_rama[rama].dtrajs)


# In[12]:

for rama in range(Nrama):
    for traj in range( len(dtraj_rama[rama]) ):
        if ( dtrajs_rama[rama][traj][np.where(dtraj_phi[rama][traj] < 0)[0][0]] != 0 ):
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == 0)[0]] -= 1
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == 1)[0]] -= 1
            dtrajs_rama[rama][traj][np.where(dtrajs_rama[rama][traj] == -1)[0]] += 2


# In[13]:

dtrajs = []
for traj in range( len(dtrajs_rama[0]) ):
    tmp = dtrajs_rama[0][traj]
    for rama in range(1,len(dtrajs_rama)):
        tmp = np.vstack( (tmp,dtrajs_rama[rama][traj]) )
    tmp = tmp.T
    dtrajs.append( tmp )
    dtrajs[traj].astype('int64')


# In[37]:

dtrajs[0].shape


# In[38]:

dtrajs_sum = []
for i in range( len(dtrajs) ):
    dtrajs_sum.append( np.sum(dtrajs[i],axis=1) )
    dtrajs_sum[i].astype('int64')


# In[39]:

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


# In[46]:

# we need a single dimensional identifier of the microstate, can we cluster to automize?
# n_clusters = 2**15
# clustering = coor.cluster_regspace(dtrajs,max_centers=n_clusters,dmin=0.5)
#clustering = coor.cluster_kmeans(dtrajs,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)

# this is unnecessary and expensive for so many states, we already generated the states, just read them in
cc = np.load('/data/isilon/rudzinski/cluster_tmp/AAQAA/AAQAA_hybrid_AMBER_Go/wDB-HP_inter/NC_CA/2016_10_21/cc_full.npy')


# In[47]:

#dtrajs_full_1D = clustering.dtrajs


# In[55]:

# cc = clustering.clustercenters[:]
cc.shape


# In[ ]:

# convert the full traj to 1D
dtrajs_full_1D = []
for traj in range(len(dtrajs)):
    dtrajs_full_1D.append([])
    for fr in range(dtrajs[traj].shape[0]):
        state = np.where( np.all(dtrajs[traj][fr]==np.array(cc),axis=1) == True )[0][0]
        dtrajs_full_1D[traj].append(state)
    dtrajs_full_1D[traj] = np.array(dtrajs_full_1D[traj])
    print 'done with traj '+str(traj)

# In[ ]:

np.save('dtrajs_rama_2st_allres_1D', dtrajs_full_1D)
#
np.save('dtrajs_Qhel_allres_1D', dtrajs_1D)


# In[ ]:

np.save('traj_rama.npy', dtraj_rama)


