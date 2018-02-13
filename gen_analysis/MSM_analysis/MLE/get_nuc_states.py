import numpy as np

cc_full = np.load('../../../../../cc_full.npy')
state_LR_weights = np.load('../../../../../state_LR_weights.npy')
mu_clust = np.load('mu_clust.npy')

Nres = cc_full.shape[1]

Nclust = mu_clust.shape[0]
state_LR_weights_clust = [[] for clust in range(Nclust)]
for clust_ind in range(Nclust):
    state_LR_weights_clust[clust_ind] = np.sum(mu_clust[clust_ind]*state_LR_weights.T,axis=1)
np.save('state_LR_weights_clust',state_LR_weights_clust)
state_LR_weights_clust = np.around(state_LR_weights_clust)

# Now, get the diff types of states on the clust level
coil = []
nuc = []
nuc2 = []
for i_ind in range(Nclust):
    if ( len(np.where(state_LR_weights_clust[i_ind]==2)[0]) == 0 ): # This is a coil state
        coil.append(i_ind)
    elif ( len(np.where(state_LR_weights_clust[i_ind]==2)[0]) == 1 ): # This is a nuc state, i.e., only a single helical segment
        nuc.append(i_ind)
    else: # check for lone helical segment res
        flag_lone_hel = False
        for res in range(1,Nres-1):
            if ( state_LR_weights_clust[i_ind][res]==2 and state_LR_weights_clust[i_ind][res-1]!=2 and state_LR_weights_clust[i_ind][res+1]!=2 ):
                flag_lone_hel = True
        if ( flag_lone_hel ):
            nuc2.append(i_ind)

hel = np.delete(np.arange(Nclust),np.hstack((coil,nuc)))

np.save('coil_states',coil)
np.save('nuc_states',nuc)
np.save('nuc2_states',nuc2)
np.save('hel_states',hel)



mu_clust_mpp = np.load('mu_clust_mpp.npy')

Nclust = mu_clust_mpp.shape[0]
state_LR_weights_clust = [[] for clust in range(Nclust)]
for clust_ind in range(Nclust):
    state_LR_weights_clust[clust_ind] = np.sum(mu_clust_mpp[clust_ind]*state_LR_weights.T,axis=1)
np.save('state_LR_weights_clust_mpp',state_LR_weights_clust)
state_LR_weights_clust = np.around(state_LR_weights_clust)

# Now, get the diff types of states on the clust level
coil = []
nuc = []
nuc2 = []
for i_ind in range(Nclust):
    if ( len(np.where(state_LR_weights_clust[i_ind]==2)[0]) == 0 ): # This is a coil state
        coil.append(i_ind)
    elif ( len(np.where(state_LR_weights_clust[i_ind]==2)[0]) == 1 ): # This is a nuc state, i.e., only a single helical segment
        nuc.append(i_ind)
    else: # check for lone helical segment res
        flag_lone_hel = False
        for res in range(1,Nres-1):
            if ( state_LR_weights_clust[i_ind][res]==2 and state_LR_weights_clust[i_ind][res-1]!=2 and state_LR_weights_clust[i_ind][res+1]!=2 ):
                flag_lone_hel = True
        if ( flag_lone_hel ):
            nuc2.append(i_ind)

hel = np.delete(np.arange(Nclust),np.hstack((coil,nuc)))

np.save('coil_states_mpp',coil)
np.save('nuc_states_mpp',nuc)
np.save('nuc2_states_mpp',nuc2)
np.save('hel_states_mpp',hel)

