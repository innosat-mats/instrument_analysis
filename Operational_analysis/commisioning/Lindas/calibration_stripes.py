#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage, calibrate_CCDitems
import datetime as DT
from selection_tools.itemselect import select_on_time as seltime
from rawdata.time_tools import add_datetime as add_datetime
from geolocation import satellite as satellite
from imagetools.imagetools import shift_image
from plotting.plotCCD import orbit_plot, simple_plot
import pandas as pd



instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/'
RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_all/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'

_,df = read_CCDdata(RacOut)
df_all=df.copy()


# %%

df=df_all


date1 = DT.datetime(2022,11,23,10,19,0)
date2 = DT.datetime(2022,11,23,10,20,0)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022



CCDSELECT=6
df = df[(df.CCDSEL == CCDSELECT)]


#%%
CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))

#%%
CCDitemsdf = pd.DataFrame.from_dict(CCDitems)
simple_plot(CCDitemsdf, image_path, nstd=2, cmap='inferno', custom_cbar=False,
    ranges=[0, 1000], format='png')



#%%

calibrate=True
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:
        totbin=CCDitem['NCBIN FPGAColumns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin

if calibrate:
    signallabel='Signal [10^10*ph/cm2/str/nm]'
    image_specification='image_calibrated'
else:
    signallabel='Counts'
    image_specification='IMAGE'

for CCDitem in CCDitems:
    #Shift image, i.e. put image on common field of view
    image_common_fov, error_flags_flipnshift = shift_image(CCDitem, CCDitem[image_specification])
    CCDitem['image_common_fov']=image_common_fov


#%%
image_specification='IMAGE'
for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    sp=plot_CCDimage(CCDitem['IMAGE'], fig, ax, title='IMAGE')
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    sp=plot_CCDimage(CCDitem['image_calibrated'], fig, ax, title='image_calibrated')
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    sp=plot_CCDimage(CCDitem['image_bias_sub'], fig, ax, title='image_bias_sub')
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    sp=plot_CCDimage(CCDitem['image_desmeared'], fig, ax, title='image_desmeared')




# %%
