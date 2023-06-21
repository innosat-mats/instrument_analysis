#%%
import matplotlib.pylab as plt
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import LSQBivariateSpline
#%%
data=loadmat('/Users/donal/Downloads/AlbedoFM_Calib_Tdep_Radiance_vs_bits.mat')

# %%
FM=data['SignFM2_Rad_raw']
temps=data['Temperatur'].squeeze()
bits=data['bitar'].squeeze()
plt.figure()
plt.pcolor(temps,bits,FM.T)
# %%
tmesh,bmesh=np.meshgrid(temps,bits)
# %%
tknots=np.linspace(min(temps),max(temps),5)
bknots=np.linspace(min(bits),max(bits),5)
spline=LSQBivariateSpline(tmesh.flatten(),bmesh.flatten(),FM.T.flatten(),tknots,bknots)    
# %%
plt.figure()
plt.pcolor(temps[::10],bits[::10],spline(temps[::10],bits[::10]).T)
# %%
