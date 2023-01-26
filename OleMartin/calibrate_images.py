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
start_time = DT.datetime(2023,1,22,17,0,0)
stop_time = DT.datetime(2023,1,22,19,5,0)

df = read_MATS_data(start_time,stop_time)
CCDitems = dataframe_to_ccd_items(df)


#%%
for i in range(len(CCDitems)):
    image_lsb, image_bias_sub, image_desmeared, image_dark_sub, image_calib_nonflipped, image_calibrated, errors = L1_calibrate(CCDitems[i],instrument)
# %%
