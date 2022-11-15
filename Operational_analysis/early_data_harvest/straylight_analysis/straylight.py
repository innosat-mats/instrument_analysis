#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate


# %%
directory='../data/Pass74/'
calibration_file='../calibration_data_linda.toml'

_,df = read_CCDdata(directory)
df = df[df.CCDSEL == 3]
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
average_image = np.zeros(CCDitems[0]["image_calibrated"].shape)
for i in range(len(CCDitems)):
    average_image = average_image + CCDitems[i]["image_calibrated"]

average_image = average_image/len(CCDitems)
average_column = np.mean(average_image[:,20:30],1)

#%%
z = np.arange(0,255)/4+50
y = 0.8e6*np.exp(-z/10)

plt.plot(y,z)
plt.plot(average_column,z)

plt.xscale('log')

# %%
