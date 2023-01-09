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
from plotting.plotCCD import orbit_plot

#from plotting.plotting_functions import collapsandplot, create_imagecube
from rawdata.time_tools import add_datetime as add_datetime
from geolocation import satellite as satellite
import pandas as pd



instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/' 
#RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_LimbPointingInNadir/'
RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_all/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'
#image_path='/Users/lindamegner/MATS/retrieval/AnalysisOutput/commissioning/LimbPointingNadir/'

_,df = read_CCDdata(RacOut)
df_all=df.copy()

df['CCDSEL'].hist()

# %%

df=df_all

#hej=df_all[(df.NRBIN==200)]

date1 = DT.datetime(2022,12,6,12,00,00)
date2 = DT.datetime(2022,12,6,14,00,00)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022

#df=df[(df.CCDSEL==1)]

df.NCSKIP.hist()
df.NRSKIP.hist()


#%%

#df=df.iloc[0::100, :]

CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))

#%%

calibrate=True
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:
        totbin=CCDitem['NCBIN FPGAColumns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin
    image_specification='image_calibrated'
else:
    image_specification='IMAGE'


for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    image=CCDitem[image_specification]
    sp=plot_CCDimage(image, fig, ax, 
    title=CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000)
    +' NRSKIP: '+str(CCDitem['NRSKIP'])
    +' NCSKIP: '+str(CCDitem['NCSKIP'])
    )




# #%%

# CCDitemsdf = pd.DataFrame.from_dict(CCDitems)

# orbit_plot(CCDitemsdf, image_path, nstd=2, cmap='inferno', custom_cbar=False,
#                ranges=[0, 1000], format='png')

# #%%
# #plt.savefig(image_path+'OperationalAnalysis2'+CCDitem['channel']+'.png', dpi=1600)

# # %%

# %%
