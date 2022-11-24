#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage

# %%
main_path='/Users/lindamegner/MATS/retrieval/git/'
directory='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/nadirtest/'
#directory=main_path+'FlightData/221122_commisioning_day1/RacOut_darkcurrent_265_266/'
calibration_file=main_path+'MATS-L1-processing/scripts/calibration_data_linda.toml'
instrument_analysis_path=main_path+'instrument_analysis/'
run_from_path=instrument_analysis_path+'Operational_analysis/commisioning/nadirtest/'
image_path=run_from_path+'images/'
_,df = read_CCDdata(directory)
df = df[(df.NCOL>500)]
#df=df.iloc[0::30, :]

#%%
CCDitems = read_CCDitems(directory,items=df.to_dict('records'))

#%%
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
#%%

calibrated=False
if calibrated:
    image_specification='image_calibrated'
else:
    image_specification='IMAGE'


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
