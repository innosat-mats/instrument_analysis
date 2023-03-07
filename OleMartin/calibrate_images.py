#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items

#%% 
calibration_file ="/home/ochri/Projects/MATS/MATS-L1-processing/scripts/calibration_data.toml"    
instrument = Instrument(calibration_file)

#%% Select on explicit time
start_time = DT.datetime(2023,2,10,19,7,20)
stop_time = DT.datetime(2023,2,10,19,7,40)

df = read_MATS_data(start_time,stop_time)
CCDitems = dataframe_to_ccd_items(df)


#%%
# for i in range(len(CCDitems)):
#     image_lsb, image_bias_sub, image_desmeared, image_dark_sub, image_calib_nonflipped, image_calibrated, errors = L1_calibrate(CCDitems[i],instrument)
# %%
#for i in range(len(CCDitems)):
I = 1
image_lsb, image_bias_sub, image_linear, image_desmeared, image_dark_sub, image_calib_nonflipped, image_calibrated, errors = L1_calibrate(CCDitems[I],instrument)

# %%
from matplotlib import pyplot as plt

plt.imshow(image_lsb,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"]+ " lsb")
plt.colorbar()
plt.show()

plt.imshow(image_bias_sub,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"] + " bias subtracted")
plt.colorbar()
plt.show()

plt.imshow(image_linear,origin='lower',clim=[0,5000])
plt.title('Channel ' + CCDitems[I]["channel"] + " linerized")
plt.colorbar()
plt.show()

plt.imshow(image_desmeared,origin='lower',clim=[-100000,100])
plt.title('Channel ' + CCDitems[I]["channel"]+ " desmeared")
plt.colorbar()
plt.show()

plt.imshow(image_dark_sub,origin='lower',clim=[-100000,100])
plt.title('Channel ' + CCDitems[I]["channel"]+ " dark subtracted")
plt.colorbar()
plt.show()

plt.imshow(image_calibrated,origin='lower',clim=[-100,10])
plt.title('Channel ' + CCDitems[I]["channel"] + " flatfied and abs calibrated")
plt.colorbar()
plt.show()

# %%
