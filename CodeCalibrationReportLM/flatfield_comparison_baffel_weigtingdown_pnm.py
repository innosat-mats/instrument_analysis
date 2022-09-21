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



dire='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/20210616_Baffletests_and_LimbHouseLightLeakage/PayloadImages/'


        

filename1=dire+'1307887925915023872_1'+'.pnm'
title1='IR1 no weight'
filename2=dire+'1307888524943252480_1'+'.pnm'
title2='IR1 60s 1000 g weight at end of baffle'

fig,ax=plt.subplots(3,1,figsize=(6,5.5))
img1 = np.float64(Image.open(filename1))  # read image
img2 = np.float64(Image.open(filename2))



    #img2_scaled=img2*img1.mean()/img2.mean()
diff=img1-img2 #_scaled
    
    
plot_CCDimage(img1, fig, ax[0],title=title1)
plot_CCDimage(img2, fig, ax[1],title=title2)
plot_CCDimage(diff, fig, ax[2],title='difference panel1-panel2')
    

fig.savefig('Flatfiewld_Weighted_Down_Baffle.jpg')