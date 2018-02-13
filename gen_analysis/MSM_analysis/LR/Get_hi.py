import numpy as np

LR_const_params = np.load('LR_const_params.npz')
Nw = LR_const_params['Nw']
Nk = LR_const_params['Nk']

hi = Nw / Nk

print 'hi = '
print hi

Nres = len(Nw)-2
fh = np.sum(Nw) / Nk
fh /= float(Nres)

print '\n'
print 'fh = '+str(fh)
