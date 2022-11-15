#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:53:36 2022

@author: Ole Martin Christensen


Code for looking at time continuity in 
measurements.
"""
#%%
from plotting.sort_images import sort_images_in_dirs, sort_images_plot
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
import matplotlib.pyplot as plt

from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage
from mats_l1_processing.instrument import Instrument
import numpy as np
import pandas as pd


orbit = 75
directory='/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/early_data_harvest/data/Pass' + str(orbit) + '/'

calibration_file='/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/early_data_harvest/calibration_data_linda.toml'

instrument = Instrument(calibration_file)

calibrate=False  

_,df = read_CCDdata(directory)

df["EXP Date"] = pd.to_datetime(df["EXP Date"])
df["TMHeaderTime"] = pd.to_datetime(df["TMHeaderTime"])




#%%
#Check how many images are lost
fig, ax = plt.subplots(1, 1)
TEXPIMS = [5.7,5.7,5.7,5.7,5.7,5.7,2]*2
TEXPIMS = [10.7,10.7,10.7,10.7,10.7,10.7,4]

import matplotlib.dates as mdates
for CCD in [1,2,3,4,5,6,7]:
    df_ccd = df.loc[df["CCDSEL"] == CCD]
    TEXPIMS_CCD = TEXPIMS[CCD-1]

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.plot(df_ccd["EXP Date"],df_ccd["EXP Nanoseconds"].diff()*1e-9-TEXPIMS_CCD,'.',label=CCD)

plt.xlabel('Time')
plt.ylabel('Diff (time)')
plt.title('Difference in EXPTIME Orbit ' + str(orbit))
ax.legend()
#plt.savefig('extime_jump_' + str(orbit) + '.png')

# fig, ax = plt.subplots(1, 1)
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
# ax.plot(df["TMHeaderTime"],df["TMHeaderNanoseconds"].diff()*1e-9,'.')
# plt.xlabel('Time')
# plt.ylabel('Diff (time)')

# fig, ax = plt.subplots(1, 1)
# df["EXP Nanoseconds"].diff().hist(bins=50)
# plt.title('EXP Nanoseconds difference')
# plt.xlabel('Timedifference')

'''
No data is lost, but there is a one second jump every now and then, this is
due to SCET beeing set with every packet. Mikael Krus has suggested 
two possible solutions (patches)

1: Fix this is post-processing
2: Set SCET time only on instrument startup (subsec is set via HW signal anyway)
'''
# %%

CCD = 1
df_ccd = df[df["CCDSEL"] == CCD]
df_ccd = df_ccd.reset_index()
TEXPIMS_CCD = TEXPIMS[CCD-1]
timediff = df_ccd["EXP Nanoseconds"].diff()*1e-9-TEXPIMS_CCD
timeshift_index = timediff[np.abs(timediff)>0.5].index
long_stab_time = np.diff(timeshift_index).argmax()

df_ccd[timeshift_index[long_stab_time-1]:timeshift_index[long_stab_time]]["EXP Nanoseconds"]


correct_time = df_ccd["EXP Nanoseconds"].loc[timeshift_index[long_stab_time]:timeshift_index[long_stab_time+1]-1]

plt.plot(np.diff(correct_time*1e-9),'.')

#%%
for i in range(0,long_stab_time):
    index_of_shift_end = timeshift_index[long_stab_time-i]
    index_of_shift_start = timeshift_index[long_stab_time-i-1]
    print(str(index_of_shift_start) + ' ' + str(index_of_shift_end))
    for j in range(index_of_shift_start,index_of_shift_end):
        print(df_ccd.at[j,"EXP Nanoseconds"])
        df_ccd.at[j,"EXP Nanoseconds"] = df_ccd.at[j,"EXP Nanoseconds"]-timediff[j]
        print(df_ccd.at[j,"EXP Nanoseconds"])

for i in range(long_stab_time+1,len(timeshift_index)-1):
    index_of_shift_end = timeshift_index[i+1]
    index_of_shift_start = timeshift_index[i]
    print(str(index_of_shift_start) + ' ' + str(index_of_shift_end))
    for j in range(index_of_shift_start,index_of_shift_end+1):
        print(df_ccd.at[j,"EXP Nanoseconds"])
        df_ccd.at[j,"EXP Nanoseconds"] = df_ccd.at[j,"EXP Nanoseconds"]-1e9
        print(df_ccd.at[j,"EXP Nanoseconds"]) 
# %%

timediff = df_ccd["EXP Nanoseconds"].diff()*1e-9-TEXPIMS_CCD
plt.plot(timediff[0:10])

# %%
