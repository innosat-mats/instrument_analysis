#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:27:10 2020

@author: lindamegner

Coompares two images from pnmfile 
"""
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

#from experimental_utils import plot_CCDimage


def plot_CCDimage(image, fig, axis, title="", clim=999):
    import numpy as np
    import matplotlib.pyplot as plt

    #sp = axis.pcolormesh(image, cmap=plt.cm.jet)
    sp = axis.imshow(image, cmap=plt.cm.jet)
    if clim == 999:
        mean = np.mean(image)
        std = np.std(image)
        sp.set_clim([mean - 1 * std, mean + 1 * std])
    else:
        sp.set_clim(clim)
        
        
    fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    return sp


# def imshowimage(fig, ax, image):
#     plt.imshow(pic)
#     if clim == 999:
#         mean = pic.mean()
#         std = pic.std()
#         plt.clim([mean - 2 * std, mean + 2 * std])
#     else:
#         plt.clim(clim)
#     plt.colorbar()

channellist=['IR1', 'IR2', 'IR3', 'IR4']

dire='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/20210616_Baffletests_and_LimbHouseLightLeakage/PayloadImages/'
for channel in channellist:

        
    #dark images images
  
    if channel=='IR1':
        filename1=dire+'1307868797249832192_1'+'.pnm'
    elif channel=='IR2':
        filename1=dire+'1307868895527298048_4'+'.pnm'
    elif channel=='IR3':
        filename1=dire+'1307868862412872192_3'+'.pnm'
    elif channel=='IR4':
        filename1=dire+'1307868830194183424_2'+'.pnm'

    
    
    
    where='upper_outer'
    
    if where=='lower_outer':
        if channel=='IR1':
            filename2=dire+'1307867722171813888_1'+'.pnm'
        elif channel=='IR2':
            filename2=dire+'1307867885441207808_4'+'.pnm'
        elif channel=='IR3':
            filename2=dire+'1307867843554763776_3'+'.pnm'
        elif channel=='IR4':
            filename2=dire+'1307867779340209920_2'+'.pnm'
    
    
    if where=='upper_inner':
        if channel=='IR1':
            filename2=dire+'1307868029010955776_1'+'.pnm'
        elif channel=='IR2':
            filename2=dire+'1307868196682968064_4'+'.pnm'
        elif channel=='IR3':
            filename2=dire+'1307868159061096192_3'+'.pnm'
        elif channel=='IR4':
            filename2=dire+'1307868119600494336_2'+'.pnm'
    

    if where=='upper_outer':
        if channel=='IR1':
            filename2=dire+'1307868302731948800_1'+'.pnm'
        elif channel=='IR2':
            filename2=dire+'1307868414600280832_4'+'.pnm'
        elif channel=='IR3':
            filename2=dire+'1307868381324691712_3'+'.pnm'
        elif channel=='IR4':
            filename2=dire+'1307868341301223680_2'+'.pnm'

    if where=='lower_inner':
        if channel=='IR1':
            filename2=dire+'1307867297839660544_1'+'.pnm'
        elif channel=='IR2':
            filename2=dire+'1307867441227081216_4'+'.pnm'
        elif channel=='IR3':
            filename2=dire+'1307867397158813440_3'+'.pnm'
        elif channel=='IR4':
            filename2=dire+'1307867349741500928_2'+'.pnm'
            
            
    if where=='side_inner':
        if channel=='IR1':
            filename2=dire+'1307868605810638336_1'+'.pnm'
        elif channel=='IR2':
            filename2=dire+'1307868718448745728_4'+'.pnm'
        elif channel=='IR3':
            filename2=dire+'1307868685566116352_3'+'.pnm'
        elif channel=='IR4':
            filename2=dire+'1307868639093826304_2'+'.pnm'


 
    fig,ax=plt.subplots(3,1,figsize=(6,5.5))
    img1 = np.float64(Image.open(filename1))  # read image
    img2 = np.float64(Image.open(filename2))


    #img2_scaled=img2*img1.mean()/img2.mean()
    diff=img1-img2 #_scaled
    
    
    plot_CCDimage(img1, fig, ax[0],title='dark')
    plot_CCDimage(img2, fig, ax[1],title='laser pointer at '+ where)
    plot_CCDimage(diff, fig, ax[2],title='difference panel1-panel2')
    
    fig.suptitle(channel)
    fig.savefig('LaserOnBaffle'+where+ channel+'.jpg')