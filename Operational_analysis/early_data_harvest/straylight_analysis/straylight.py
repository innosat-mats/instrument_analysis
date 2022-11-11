#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np


# %%
directory='data/Pass74/'
_,df = read_CCDdata(directory)
df = df[df.CCDSEL == 1]
CCDitems = read_CCDitems(directory,items=df.to_dict('records'))

#%%
average_image = np.zeros(CCDitems[0]["IMAGE"].shape)
for i in range(len(CCDitems)):
    average_column = average_column + CCDitems[i]["IMAGE"]

average_image = average_column/len(CCDitems)
plt.plot(average_image[:,1]-200,np.arange(0,len(average_image[:,1])))
plt.xscale('log')

# %%
plt.pcolor(average_image)
# %%
