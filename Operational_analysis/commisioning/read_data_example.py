# %%
from mats_l1_processing.experimental_utils import plot_CCDimage
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
import rawdata.time_tools as time_tools
import geolocation.satellite as satellite
import selection_tools.itemselect as itemselect
import datetime as DT

# %%
directory='/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/commisioning/data/star/'
calibration_file='/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/commisioning/calibration_data_linda.toml'

_,df = read_CCDdata(directory)

#%%
df = time_tools.add_datetime(df) #convert to datetime
start_date = DT.datetime(2022,12,5,18,00,00)
end_date = DT.datetime(2022,12,5,20,00,00)
df = itemselect.select_on_time(start_date,end_date,df)




#%%
#Read in data

CCDitems = read_CCDitems(directory,items=df.to_dict('records'))

#%% 
#Calibrate data

instrument = Instrument(calibration_file)

for CCDitem in CCDitems:
    (
        image_lsb,
        image_bias_sub,
        image_desmeared,
        image_dark_sub,
        image_calib_nonflipped,
        image_calibrated,
        errors
    ) = L1_calibrate(CCDitem, instrument)
#%% Do some processing
average_image_calib = np.zeros(CCDitems[0]["image_calibrated"].shape)
average_image = np.zeros(CCDitems[0]["IMAGE"].shape)

for i in range(len(CCDitems)):
    average_image_calib = average_image_calib + CCDitems[i]["image_calibrated"]
    average_image = average_image + CCDitems[i]["IMAGE"]

average_image_calib = average_image_calib/len(CCDitems)
average_column_calib = np.mean(average_image_calib[:,20:21],1)
average_image = average_image/len(CCDitems)
average_column = np.mean(average_image[:,20:21],1)
#%% Plot results
z = np.arange(0,255)/4+50
y = average_column_calib[80]*np.exp(-(z-70)/8.5)

plt.plot(y,z)
plt.plot(average_column_calib,z)
plt.plot(average_column,z)

plt.legend(['rayleigh','calib','raw'])
plt.xscale('log')
plt.show()
# %%
#Plot a single image
fig , ax= plt.subplots(1, 1, frameon=False, figsize=(4, 1))
sp=plot_CCDimage(CCDitems[0]['image_calibrated'], fig, ax)
plt.show()
