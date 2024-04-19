import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(0,0.1,20)
r_freeSurf = np.linspace(0.05,0.1,20)
y = np.linspace(0,0.1,20)
r_crys = np.linspace(0,0.05,20)

# Generate data
data_side = np.random.rand(y.shape[0],10) * 2000 - 1000
data_side[:,0] = 10
data_side[:,1] = 1
data_side[:,2] = np.linspace(1,20,20)
data_side[:,3] = r.max()
data_side[:,4] = y
data_side[:,5] = 0
p_gradTr = [1.90476190e+05, -4.35714286e+04,  1.80238095e+03, -5.00000000e+00]
data_side[:,6] = 1685 + (0.1-y)*50 # temperature
data_side[:,7] = np.polyval(p_gradTr, y)/2 # gradient r
data_side[:,8] = 0 # gradient y, not used
np.savetxt('side.dat', data_side, fmt='%f', delimiter=' ')

# Generate data
data_freeSurf = np.random.rand(y.shape[0],10) * 2000 - 1000
data_freeSurf[:,0] = 10
data_freeSurf[:,1] = 1
data_freeSurf[:,2] = np.linspace(1,20,20)
data_freeSurf[:,3] = r_freeSurf
data_freeSurf[:,4] = y.max()
data_freeSurf[:,5] = 0
data_freeSurf[:,6] = 1710 + r*100 # temperature
data_freeSurf[:,7] = 0  # gradient r (not used)
p = [54769.35105551, -8316.06724003,   269.571931]
data_freeSurf[:,8] =  np.polyval(p_gradTr, r_freeSurf) # gradient y
np.savetxt('freeSurf.dat', data_freeSurf, fmt='%f', delimiter=' ')

# Generate data
data_crysSurf = np.random.rand(y.shape[0],10) * 2000 - 1000
data_crysSurf[:,0] = 10
data_crysSurf[:,1] = 1
data_crysSurf[:,2] = np.linspace(1,20,20)
data_crysSurf[:,3] = r_crys
data_crysSurf[:,4] = y.max()
data_crysSurf[:,5] = 0
data_crysSurf[:,6] = 1685 # temperature
data_crysSurf[:,7] = 0 # gradient r
data_crysSurf[:,8] = 0 # gradient y
print(data_crysSurf)
print(data_crysSurf.shape)
np.savetxt('crysSurf.dat', data_crysSurf, fmt='%f', delimiter=' ')
