#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
import numpy as np
import matplotlib.pylab as plt
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
import io
from tangentlib import *
#%matplotlib widget

#%% 
calibration_file ="/Users/donal/projekt/SIW/instrument_analysis/Operational_analysis/commisioning/Donal/data/calibration_data.toml"    
instrument = Instrument(calibration_file)

#filename='/Users/donal/projekt/SIW/MATS-L1-processing/calibration_data/pointing/qprime.csv'
#qprimes=np.loadtxt(filename,delimiter=',',dtype={'names':('Channel','q0','q1','q2','q3'),
#                           'formats':('S3','f4','f4','f4','f4')})
#%%
def heights(ccditem,x):
    d = ccditem['EXP Date']
    ts = sfapi.load.timescale()
    t = ts.from_datetime(d)
    ecipos = ccditem['afsGnssStateJ2000'][0:3]
    q = ccditem['afsAttitudeState']
    quat = R.from_quat(np.roll(q, -1))
    qprime=R.from_quat(ccditem['qprime'])
    ypixels = np.linspace(0, ccditem['NROW'], 40)
    ths=np.zeros_like(ypixels)
    xdeg,ydeg=pix_deg(ccditem,x,ypixels)
    for iy,y in enumerate(ydeg):
        los = R.from_euler('XYZ', [0, y, xdeg], degrees=True).apply([1, 0, 0])
        ecivec = quat.apply(qprime.apply(los))
        ths[iy] = findtangent(t, ecipos, ecivec).fun
    return CubicSpline(ypixels,ths/1000)
# %%
def pix_deg(ccditem, xpixel, ypixel):
    """ 
    Inputs : CCDitem  , xpixel (can be an array), ypixel (can be an array)
    Outputs xdeg - x offset of the pixel(s) in degrees, ydeg - y offset of the pixels(s) in degrees
    
    """
    xdisp = 6.06/2047
    ydisp = 1.52/510  # 1.52/510
    ncskip = ccditem['NCSKIP']
    ncbin = ccditem['NCBIN CCDColumns']
    nrskip = ccditem['NRSKIP']
    nrbin = ccditem['NRBIN']
    ncol = ccditem['NCOL']
    # flipped configuration: clips from right before flip left clip by limiting no. of  columns
    if (ccditem['CCDSEL']) in [1, 3, 5, 6]:
        xdeg = xdisp*((2048-ncskip - (ncol+1)*ncbin +
                      ncbin*(xpixel+0.5)) - 2047./2)
    else:
        xdeg = xdisp*(ncskip + ncbin * (xpixel+0.5) - 2047./2)
    ydeg = ydisp*(nrskip + nrbin * (ypixel+0.5) - 510./2)
    return xdeg, ydeg

    
#%% Select on explicit time
start_time = DT.datetime(2023,1,27,17,0,0)
stop_time = DT.datetime(2023,1,27,20,5,0)
start_time = DT.datetime(2023,1,28,6,0,0)
stop_time = DT.datetime(2023,1,28,11,5,0)
start_time = DT.datetime(2023,2,6,0,0,0)
stop_time = DT.datetime(2023,2,6,11,5,0)
df = read_MATS_data(start_time,stop_time,version='0.4')
CCDitems = dataframe_to_ccd_items(df)
#%%
ccdnames=('IR1','IR4','IR3','IR2','UV1','UV2','NADIR')
ir1=dataframe_to_ccd_items(df[df.CCDSEL==1])
ir4=dataframe_to_ccd_items(df[df.CCDSEL==2])
ir3=dataframe_to_ccd_items(df[df.CCDSEL==3])
ir2=dataframe_to_ccd_items(df[df.CCDSEL==4])
uv1=dataframe_to_ccd_items(df[df.CCDSEL==5])
uv2=dataframe_to_ccd_items(df[df.CCDSEL==6])
#%%
channels=[ir1,ir2,ir3,ir4,uv1,uv2]
#for ch in channels:
#    qprime=[(d['q0'],d['q1'],d['q2'],d['q3']) for d in qprimes if d['Channel'].decode('utf8')==ccdnames[ch[0]['CCDSEL']-1]]
#    for ccditem in ch:
#        ccditem['qprime']=np.array(qprime[0])
#%%
plt.figure(figsize=[4,2])
#image=ir1[n]['IMAGE']
image=ir1[0]['IMAGE']
plt.imshow(image,origin='lower')
[col, row]=image.shape
mean = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].mean()
std = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].std()
plt.clim([mean - 2 * std, mean + 2 * std])
plt.gca().axis("auto")
plt.colorbar()
#%%
def calibrate(CCDitem, instrument):
    (image_lsb, image_bias_sub, 
    image_desmeared, image_dark_sub, 
    image_calib_nonflipped, image_calibrated, errors) = L1_calibrate(CCDitem,instrument)
    return image_calibrated
#%%
#n=494 #wird nlc
n=0
n=94
plt.close('all')
ir1cal=calibrate(ir1[n],instrument)
ir2cal=calibrate(ir2[n],instrument)
ir3cal=calibrate(ir3[n],instrument)
ir4cal=calibrate(ir4[n],instrument)
uv1cal=calibrate(uv1[n],instrument)
uv2cal=calibrate(uv2[n],instrument)
# %%
plt.figure(figsize=[4,2])
plt.imshow(uv1cal,origin='lower')
plt.gca().axis("auto")
plt.clim([0,600])
plt.colorbar()
# %%
plt.figure(figsize=[4,2])
plt.semilogy(ir1cal[:,24]-0)
plt.semilogy(ir2cal[:,24]-0)
plt.figure(figsize=[4,2])
plt.semilogy(ir3cal[:,4]-0)
plt.semilogy(ir4cal[:,4]-0)
  # %%
print(ir1[n]['EXP Date'])
print(ir2[n]['EXP Date'])


# %%
csh=heights(ir4[n],8)
# %%
plt.figure(figsize=[6,4])
csh=heights(ir1[n],24+6)
y4=csh(range(ir1cal.shape[0]))
plt.semilogx(ir1cal[:,24+6]-0,y4)
csh=heights(ir2[n],24)
y2=csh(range(ir2cal.shape[0]))
plt.semilogx(ir2cal[:,24]-0,y2)
plt.ylabel('Altitude(km)')
plt.xlabel('Intensity (Witt)')
plt.legend(('IR1','IR2'))

plt.figure(figsize=[6,4])
csh=heights(ir3[n],4)
y3=csh(range(ir3cal.shape[0]))
plt.semilogx(ir3cal[:,4]-0,y3)
csh=heights(ir4[n],4)
y4=csh(range(ir4cal.shape[0]))
plt.semilogx(ir4cal[:,4]-0,y4)
plt.gca().axis('auto')
plt.autoscale()
plt.ylabel('Altitude(km)')
plt.xlabel('Intensity (Witt)')
plt.legend(('IR3','IR4'))


# %%
# %%
print(ir1[n]['EXP Date'])
print(ir2[n]['EXP Date'])
print(ir3[n]['EXP Date'])
print(ir4[n]['EXP Date'])
# %%
uv1[1]['qprime']
# %%
uv1cal_0=calibrate(uv1[0],instrument)
uv1cal_1=calibrate(uv1[1],instrument)
plt.figure()
data=uv1cal_0[323-150:323+30,1044-44:1044+44]-uv1cal_1[245-150:245+30,1299-44:1299+44]
plt.imshow(data,origin='lower')
plt.clim([-200,200])
plt.colorbar()

# %%
plt.figure(figsize=[8,2])
plt.imshow(uv1cal_0,origin='lower')
plt.gca().axis("auto")
plt.clim([0,5000])
plt.colorbar()
# %%
