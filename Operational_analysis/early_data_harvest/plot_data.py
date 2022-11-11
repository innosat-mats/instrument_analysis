#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems
from matplotlib import pyplot as plt

directory='data/Pass70/'
CCDitems = read_CCDitems(directory)


fig, ax = plt.subplots(1)
plotCCDitem(CCDitems[0], fig, ax, title='First Image !!!')
ax.set_aspect('equal')
# %%
