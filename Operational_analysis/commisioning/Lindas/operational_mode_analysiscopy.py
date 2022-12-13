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

# %%

instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/'
RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'
# %%
_,df = read_CCDdata(RacOut)
df.columns = [c.replace(' ', '_') for c in df.columns]
df_all=df
# %%
date1 = DT.datetime(2022,11,23,12,00,00)
date2 = DT.datetime(2022,11,25,12,00,00)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022

#%%

df = df[(df.CCDSEL == 2) & (df.NCOL < 500) ]
df=df.iloc[0:10, :]

#%%
CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))


#%%
calibrate=True
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    image_specification='image_calibrated'
else:
    image_specification='IMAGE'



#%%

clim=[-5, 100]
meanvalue=[]
exptime=[]
for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    image=CCDitem[image_specification]
    sp=plot_CCDimage(image, fig, ax, title=CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000))
    meanvalue.append(image.mean())
    exptime.append(CCDitem['TEXPMS']/1000)


#%%
    
plt.plot(exptime,meanvalue,'*')
plt.xlabel('Integration time')
plt.ylabel('Mean signal')
plt.title(CCDitem['channel'])
plt.show()




#plt.savefig(image_path+CCDitem['id']+'.png', dpi=1600)

# %%
