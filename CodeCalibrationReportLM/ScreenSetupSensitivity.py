#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""





from mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage,read_all_files_in_protocol 
from mats_l1_processing.read_in_functions import readprotocol 
import matplotlib.pyplot as plt





#directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/Diffusor/DiffusorFlatTests/'
#protocol='ABOUT.txt'

# directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200330_flatfields_0C/'
# protocol='flatfields_200330_SigMod0_LMprotocol.txt'

directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/201005_ScreenSensitivity/'
protocol='protocol.txt'

read_from='rac'  
df_protocol=readprotocol(directory+protocol)

#df_bright=df_protocol[df_protocol.DarkBright=='B']
CCDItems=read_all_files_in_protocol(df_protocol, read_from,directory)



channel='IR1'
# if channel=='IR1':
#      longtime=14000
# elif channel == 'UV2':
#      longtime=60000
# shorttime=2000     
     
Brights= list(filter(lambda x: (x['DarkBright']=='B'),CCDItems))
Darks= list(filter(lambda x: ( x['DarkBright']=='D'),CCDItems))

pictures={}
for i, CCDitem in enumerate(Darks):

    pictures[i]=Brights[i]['IMAGE']-Darks[i]['IMAGE']
    
    
fig, ax =plt.subplots(3,2,figsize=(10,8))

    #plotCCDitem(CCDitem,fig, ax,CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000)+'s')   

plot_CCDimage(pictures[0], fig, ax[0,0], 'original 28.9')
plot_CCDimage(pictures[1], fig, ax[1,0], '27.7')
plot_CCDimage(pictures[2], fig, ax[2,0], 'back to original 28.9')


# meanpicture=(pictures[0]+pictures[1]+pictures[2])/3.

# plot_CCDimage((pictures[0]-meanpicture), fig, ax[0,1], 'diff mean')
# plot_CCDimage((pictures[1]-meanpicture), fig, ax[1,1], 'diff mean')
# plot_CCDimage((pictures[2]-meanpicture), fig, ax[2,1], 'diff mean')


plot_CCDimage(pictures[0]-pictures[0], fig, ax[0,1], 'diff original')
plot_CCDimage(pictures[1]-pictures[0], fig, ax[1,1], 'diff original')
plot_CCDimage(pictures[2]-pictures[0], fig, ax[2,1], 'diff original')


fig.savefig('ScreenSetupSensitivity'+channel+'.jpg')
