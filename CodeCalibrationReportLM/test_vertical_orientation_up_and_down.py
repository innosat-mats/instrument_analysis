#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 14:04:08 2022

@author: lindamegner
"""



from mats_l1_processing.experimental_utils import  plotCCDitem, read_all_files_in_protocol, plot_CCDimage
from mats_l1_processing.experimental_utils import readprotocol

from mats_l1_processing.L1_calibration_functions import (
    get_true_image,
    desmear_true_image,
    CCD,
    subtract_dark)

import matplotlib.pyplot as plt



directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/LaserAlignmentTest_RacFiles0909to0910Duplicated/'
directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/LegoTests/LegoTests210913/'

protocol='made_up_protocol_to_test_orientation.txt'


read_from="rac" 


fig, ax=plt.subplots(4,1)
df_protocol=readprotocol(directory+protocol)
CCDitems_all=read_all_files_in_protocol(df_protocol, read_from,directory)

channel='UV2'
CCDitems=[]
for CCDitem in CCDitems_all:
    if CCDitem['channel']==channel:
        CCDitems.append(CCDitem)


plotCCDitem(CCDitems[0],fig, ax[0], title=CCDitems[0]['channel']+ ' '+CCDitems[0]['id']+'dark')
plotCCDitem(CCDitems[1],fig, ax[1], title=CCDitems[1]['channel']+ ' '+CCDitems[1]['id']+'fully lit')
plotCCDitem(CCDitems[2],fig, ax[2], title=CCDitems[2]['channel']+ ' '+CCDitems[2]['id']+'lego inserted at the lower edge')
plot_CCDimage(CCDitems[2]['IMAGE']-CCDitems[1]['IMAGE'], fig, ax[3],title='panel 2 - panel 3' )



#for i, CCDitem in enumerate(CCDitems[0::5]):    
#    plotCCDitem(CCDitem,fig, ax[i], title=CCDitem['channel']+ ' '+CCDitem['id'])
     


fig.savefig('images/vertica_orientation.jpg') 
    

CCDitem=CCDitems[0]
image=CCDitem['IMAGE']
image_bias_sub,error_flags_bias = get_true_image(CCDitem)


fig, ax=plt.subplots(2,1)
plotCCDitem(CCDitem,fig, ax[0], title=CCDitem['channel']+ ' '+'dark')
plot_CCDimage(image_bias_sub,fig, ax[1], title=CCDitem['channel']+ ' '+'bias subtracted')

exptime=20000
mysignal=image_bias_sub.mean()
mean_bias=image.mean()-mysignal
expected_value=mysignal*exptime/CCDitem['TEXPMS']+mean_bias
print('expected value at '+str(exptime)+'s = '+ str(expected_value))

print('mean bias:  '+str(mean_bias))


