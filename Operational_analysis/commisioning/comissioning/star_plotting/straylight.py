#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate


# %%
directory='/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/comissioning/star_plotting/star_OHB/'
calibration_file='../calibration_data_linda.toml'

_,df = read_CCDdata(directory)
df = df[df.CCDSEL == 3]
CCDitems = read_CCDitems(directory,items=df.to_dict('records'))

plt.imshow(CCDitems[-1]['IMAGE'])
plt.show()