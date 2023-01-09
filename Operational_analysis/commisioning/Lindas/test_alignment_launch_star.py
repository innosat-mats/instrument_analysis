#%% 
 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 21:00:16 2022

@author: lindamegner


"""



from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import calibrate_CCDitems
import datetime as DT
from selection_tools.itemselect import select_on_time as seltime
from rawdata.time_tools import add_datetime as add_datetime
from geolocation import satellite as satellite
from imagetools.imagetools import shift_image
from plotting.plotCCD import orbit_plot, simple_plot
import pandas as pd

import warnings
#from mats_l1_processing.L1b_calibration_functions import shift_image

def plot_CCDimage(image, fig, axis, title="", clim=999, aspect="auto"):
    sp = axis.imshow(image, cmap="magma", origin="lower", interpolation="none")
    # sp=axis.pcolormesh(image, , cmap='viridis')
    if clim == 999:
        [col, row]=image.shape
        #Take the mean and std of the middle of the image, not boarders
        mean = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].mean()
        std = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].std()
        sp.set_clim([mean - 2 * std, mean + 2 * std])
    else:
        sp.set_clim(clim)
    fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    axis.set_aspect(aspect)
    return sp

def plot_CCDimage_transp(image, fig, axis, title="", clim=999, aspect="auto", alpha=1, 
    ncshift=0, nrshift=0, width=0, height=0, pos00=[0,0], colour='black'):
    sp = axis.imshow(image, cmap="viridis", origin="lower", interpolation="none", alpha=alpha)
    # sp=axis.pcolormesh(image, , cmap='viridis')
    if clim == 999:
        mean = image.mean()
        std = image.std()
        sp.set_clim([mean - 2 * std, mean + 2 * std])
    else:
        sp.set_clim(clim)
    #fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    rectangle = plt.Rectangle((ncshift, nrshift), width, height, facecolor='none', ec=colour)
    axis.add_patch(rectangle)
    #rectangle = plt.Rectangle((pos00[1], pos00[0]), 2047, 511, facecolor='none', ec=colour)
    #axis.add_patch(rectangle)
    axis.set_aspect(aspect)
    
    return sp
#%%
def get_crop_positions(pos00, channel):
    class channelinfo:
       def __init__(self, x_pos, y_pos):
            self.x_pos=x_pos
            self.y_pos=y_pos

    info_IR1=channelinfo(x_pos=-75,y_pos=47)
    info_IR2=channelinfo(x_pos=144,y_pos=76)
    info_IR3=channelinfo(x_pos=37,y_pos=66)
    info_IR4=channelinfo(x_pos=0,y_pos=0)
    info_UV1=channelinfo(x_pos=88,y_pos=13)
    info_UV2=channelinfo(x_pos=156,y_pos=192)

    nrowmax=511
    ncolmax=2047
    ncol=ncolmax
    nrow=nrowmax

    ncskip=0
    nrskip=0
    if channel=='IR1':
        ncskip=info_IR2.x_pos-info_IR1.x_pos
        ncol=ncolmax-(info_IR2.x_pos-info_IR1.x_pos)
        nrskip=0
        nrow=nrowmax-(info_IR2.y_pos-info_IR1.y_pos)
    elif channel=='IR2':
        ncskip=0
        ncol=ncolmax-(info_IR2.x_pos-info_IR1.x_pos)
        nrskip=0
        nrow=nrowmax   
    elif channel=='IR3':
        ncskip=info_IR2.x_pos-info_IR3.x_pos
        ncol=ncolmax-(info_IR2.x_pos-info_IR1.x_pos)
        nrskip=info_IR3.y_pos-info_IR4.y_pos
        nrow=nrowmax-(info_IR2.y_pos-info_IR3.y_pos)
    elif channel=='IR4':
        ncskip=info_IR2.x_pos-info_IR4.x_pos
        ncol=ncolmax-(info_IR2.x_pos-info_IR1.x_pos)
        nrskip=0
        nrow=nrowmax-(info_IR2.y_pos-info_IR4.y_pos)
    elif channel=='UV1':
        ncskip=info_UV2.x_pos-info_UV1.x_pos
        ncol=ncolmax-info_UV2.x_pos+info_UV1.x_pos
        nrskip=0
        nrow=nrowmax-(info_UV2.y_pos-info_UV1.y_pos)
    elif channel=='UV2':
        ncskip=0
        ncol=ncolmax-info_UV2.x_pos+info_UV1.x_pos
        nrskip=info_UV2.y_pos-info_UV1.y_pos
        nrow=nrowmax-(info_UV2.y_pos-info_UV1.y_pos)



    ncshift=pos00[1]+ncskip
    nrshift=pos00[0]+nrskip

    width=ncol
    height=nrow
    return ncshift,nrshift,ncskip,nrskip, width, height


star=999

instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/'
if star==999: 
    RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_November/'
else:
    RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_all/'
#RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_from25nov/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'


#%%
_,df = read_CCDdata(RacOut)
df_all=df.copy()


# %%

df=df_all


if star==4:
    date1 = DT.datetime(2022,11,29,8,30,00)
    date2 = DT.datetime(2022,11,29,10,30,00)
elif star==5:
    date1 = DT.datetime(2022,11,29,10,30,00)
    date2 = DT.datetime(2022,11,29,12,00,00)
elif star==1:
    date1 = DT.datetime(2022,11,29,12,00,00)
    date2 = DT.datetime(2022,11,29,13,00,00)
elif star==6:
    date1 = DT.datetime(2022,11,29,13,00,00)
    date2 = DT.datetime(2022,11,29,15,00,00)
elif star==7000:
    date1 = DT.datetime(2022,12,12,21,00,00)
    date2 = DT.datetime(2022,12,12,23,00,00)
elif star==7:
    date1 = DT.datetime(2022,12,5,19,50,00)
    date2 = DT.datetime(2022,12,5,20,30,00)
elif star==999: #full frame staring mode test LBFT
    date1 = DT.datetime(2022,11,22,14,19,00)
    date2 = DT.datetime(2022,11,23,14,51,00)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022


CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))

# #%%
# CCDitemsdf = pd.DataFrame.from_dict(CCDitems)
# simple_plot(CCDitemsdf, image_path, nstd=2, cmap='inferno', custom_cbar=False,
#     ranges=[0, 1000], format='png')

#%%
image_specification='IMAGE'
for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    image=CCDitem[image_specification]
    sp=plot_CCDimage(image, fig, ax, title=CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000))


#%%



calibrate=True
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:        
        totbin=CCDitem['NCBIN FPGAColumns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin
        #Shift image, i.e. put image on common field of view
        image_common_fov, error_flags_flipnshift = shift_image(CCDitem, CCDitem['image_calibrated'])
        CCDitem['image_common_fov']=image_common_fov

if calibrate:
    signallabel='Signal [10^10*ph/cm2/str/nm]'
    image_specification='image_calibrated'
else:
    signallabel='Counts'
    image_specification='IMAGE'



#%%

fig1, ax1 = plt.subplots(6, 1)
fig2, ax2 = plt.subplots(1)
fig3, ax3 = plt.subplots(6, 1)
start=0
stop=5
fix_clim=[-5,20]
colours=['orange', 'blue', 'red', 'yellow','cyan','magenta', 'white']
for ind, CCDitem in enumerate(CCDitems):
    colour=colours[ind]
    plot_CCDimage(CCDitem['image_calibrated'], fig1, ax1[ind], CCDitem['channel'])    
    plot_CCDimage(CCDitem['image_common_fov'], fig3, ax3[ind], title=CCDitem['channel'], clim=fix_clim)  

    pos00=np.argwhere(np.isfinite(CCDitem['image_common_fov']))[0]
    ncshift, nrshift, ncskip,nrskip,width, height=get_crop_positions(pos00, CCDitem['channel'])
    CCDitem['set_ncskip']=ncskip
    CCDitem['set_nrskip']=nrskip
    CCDitem['set_ncol']=width
    CCDitem['set_nrow']=height

    sp=plot_CCDimage_transp(CCDitem['image_common_fov'], fig2, ax2, 
    title='All channels', clim=fix_clim, alpha=0.3, ncshift=ncshift, nrshift=nrshift, 
        width=width, height=height, pos00=pos00, colour=colour)    
    
    ax2.text(pos00[1], pos00[0], str(CCDitem['channel']) , color=colour)
    rectangle = plt.Rectangle((ncshift, nrshift), width, height, facecolor='none', ec=colour)
    #ax2.add_patch(rectangle)
    ax3[ind].add_patch(rectangle)
    

    
    fig2.savefig('images/Alignment_crop_together_'+str(star)+'.png',  dpi=700) 
    fig3.savefig('images/Alignment_crop_all_'+str(star)+'.png',  dpi=700)   
  

for CCDitem in CCDitems:
    if CCDitem['channel']=='IR1' or CCDitem['channel']=='IR3' or CCDitem['channel']=='UV1' or CCDitem['channel']=='UV2':
        CCDitem['set_ncskip']=2047-(CCDitem['set_ncol']+CCDitem['set_ncskip'])
    print (CCDitem['channel']+' NCSKIP: '+str(CCDitem['set_ncskip'])
        +' NRSKIP: '+str(CCDitem['set_nrskip'])
        +' NCOL: '+str(CCDitem['set_ncol'])
        +' NROW: '+str(CCDitem['set_nrow']))



#fig.savefig('Testfile.jpg')    

#  rectangle = plt.Rectangle((10,10), 200, 200, fc='blue',ec="red")
  #  plt.gca().add_patch(rectangle) 

# %%
