#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 10:13:48 2021

@author: lindamegner
"""

from mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage, read_all_files_in_protocol, plotCCDitem
from mats_l1_processing.read_in_functions import readprotocol
import matplotlib.pyplot as plt
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.L1_calibration_functions import (
    get_true_image,
    get_linearized_image,
    desmear_true_image,
    CCD,
    subtract_dark,
    compensate_flatfield,
    get_linearized_image,
)


def subtract_images(CCDitem1, CCDitem2):
    diffimage=(CCDitem1['IMAGE']-CCDitem2['IMAGE'])
    return diffimage

    


directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Polarisation_and_Relative_Calibration/'
protocol='filelist_polarisation_IR_no_comma.txt'
#protocol='filelist_polarisation_UV_no_comma.txt'


read_from="rac" 
df_protocol=readprotocol(directory+protocol)


CCDitems=read_all_files_in_protocol(df_protocol[:10], read_from,directory)






#Example of how to plot (3 first images)
for CCDitem in CCDitems[:7]:
   fig=plt.figure()
   ax=fig.gca()
   plotCCDitem(CCDitem,fig, ax,CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000)+'s') 
   print(CCDitem['id'])
   
   
#Example of how to plot using only the image, not CCDitems:
fig1=plt.figure()
ax1=plt.subplot(2,1,1)
diffimage=subtract_images(CCDitems[0],CCDitems[1])
plot_CCDimage(diffimage, fig1, ax1, title='hej jörg')
#OR 


ax2=plt.subplot(2,1,2)
sp=ax2.imshow(diffimage, cmap="viridis", origin="lower", interpolation="none")
fig1.colorbar(sp, ax=ax2)
ax2.set_title('second way to plot diffference')
ax2.set_aspect("auto")




# Now desmear
nr_of_pics=5
fig,ax=plt.subplots(nr_of_pics,1)
for i, CCDitem in enumerate(CCDitems[:nr_of_pics]):

    #image_linear = get_linearized_image(CCDitem, diffimage)

     image_desmeared = desmear_true_image(CCDitem, CCDitem['IMAGE'])


     plot_CCDimage(image_desmeared,fig, ax[i], 'Desmeared')  