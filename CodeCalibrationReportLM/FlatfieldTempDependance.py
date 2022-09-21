#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:25:11 2020

@author: lindamegner
"""


from mats_l1_processing.items_units_functions import read_files_in_protocol_as_ItemsUnits 
from mats_l1_processing.experimental_utils import plot_CCDimage_hmean
from mats_l1_processing.experimental_utils import readprotocol
import matplotlib.pyplot as plt
import numpy as np


directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200511_temperature_dependence/'



read_from='rac'  

df_protocol=readprotocol(directory+'protocol.txt')
#df_only2 = df_protocol[(df_protocol.index-2) % 3 != 0]

CCDItemsUnits=read_files_in_protocol_as_ItemsUnits(df_protocol,directory,3,read_from)

fig_all=plt.figure()
ax_all=fig_all.gca()


myfig, myax=plt.subplots()
fig_brus, ax_brus=plt.subplots()
mybrusfig, mybrusax=plt.subplots()



channels=['IR1','IR2','IR3','IR4','UV1','UV2']
for ichan, channel in enumerate(channels[:]):
    CCDItemsUnitsSelect= list(filter(lambda x: ( x.imageItem['channel']==channel),CCDItemsUnits))
    
    
    if channel[0]=='I':
        maxplot=8
    elif channel[0]=='U':
        maxplot=7
    else:
        raise Exception('Unknown channel name')
    fig,ax=plt.subplots(maxplot,1)
    imgmean=[]
    imgmeanerr=[]
    temperature=[]
    brus=[]
    for i, ItemsUnit in enumerate(CCDItemsUnitsSelect[0:]):
         # fig=plt.figure()
         # ax=fig.gca()    
         ItemsUnit.plot(fig,ax[i], title=ItemsUnit.imageItem['channel'])
         ax[i].text(100,200,'mean:'+str(np.mean(ItemsUnit.subpic)))
         ax[i].text(100,320,'T:'+str(ItemsUnit.imageItem['temperature']))
         subpixmid=ItemsUnit.subpic[100:400,100:1948]
         nr_of_elements=subpixmid.size
         
         if True: # Mask outsiders
             low_y = (np.mean(subpixmid) - np.std(subpixmid))
             high_y=(np.mean(subpixmid) + np.std(subpixmid))         
             subpixmid_mask = np.ma.masked_outside(subpixmid,low_y,high_y)
             subpixmid=subpixmid_mask
             nr_of_elements=subpixmid_mask.count()
         

         imgmean.append(np.mean(subpixmid))
         imgmeanerr.append(np.std(subpixmid)/np.sqrt(nr_of_elements))
         temperature.append(ItemsUnit.imageItem['temperature'])
         brus.append(np.std(ItemsUnit.dark1Item['IMAGE']-ItemsUnit.dark2Item['IMAGE']))
    fig.suptitle('Flat fields as function of temperature')
   # fig.savefig('FlatfeildTempVariability_200617_'+channel+'.jpg')
    

    imgmeanmean=np.mean(imgmean)
    
    
    ax_all.errorbar(temperature, imgmean/imgmean[0],yerr=imgmeanerr/imgmean[0], label=channel)
    
    
    
    #myplot=myax.plot(temperature[:], imgmean[:]/imgmean[0], marker='*', label=channel)
    #mycolor = myplot[0].get_color()

    z, res, _, _, _ = np.polyfit(temperature[:], imgmean[:]/imgmean[0], 1, full=True)
    p = np.poly1d(z)
    xp = np.linspace(5, 23, 100)
    reftemp=13

    if channel=='IR3':
        channelinfo=channel+ ' (old)'
    else:
        channelinfo=channel
        
        
    
    
    myplot=myax.plot(temperature[:], imgmean[:]/imgmean[0]/p(reftemp), marker='*', label=channelinfo)
    mycolor = myplot[0].get_color()
    
    myax.plot(xp, p(xp)/p(reftemp), '-',  color=mycolor)
    myax.text(9,1.025-0.003*ichan, 'k='+str.format('{0:.7f}', z[0])+' /$^\circ$C', color=mycolor)
    #myax.text(12,1.02-0.003*ichan, 'res='+str.format('{0:.4g}', res[0]), color=mycolor)
    plot_brus=ax_brus.plot(temperature, brus/brus[0], marker='*', label=channelinfo, color=mycolor)

    #Now plot as a function of noise level instead
    z, res, _, _, _ = np.polyfit(brus[:]/np.mean(brus), imgmean[:]/np.mean(imgmean), 1, full=True)
    p = np.poly1d(z)
    xp = np.linspace(0.8, 1.3, 100)
    reftemp=6
    mybrusplot=mybrusax.plot(brus[:]/np.mean(brus), imgmean[:]/np.mean(imgmean), marker='*', color=mycolor, label=channelinfo)
    mybrusax.plot(xp, p(xp), '-',  color=mycolor)
    mybrusax.text(.8,1.06-0.01*ichan, 'k='+str.format('{0:.3}', z[0]), color=mycolor)


myax.set_xlabel('Temperature [$^\circ$C]')
myax.set_ylabel('Mean signal strength relative to ' +str(reftemp)+'$^\circ$C')
myfig.show()  
myax.legend()  
myfig.savefig('signal_strength_vs_temperature.jpg')


ax_brus.set_xlabel('Temperature [$^\circ$C]')
ax_brus.set_ylabel('Relative noise level')
ax_brus.legend()  
fig_brus.savefig('temperature_vs_noise.jpg')

mybrusax.set_xlabel('Noise [LSB]')
mybrusax.set_ylabel('Mean relative signal strength') 
mybrusax.legend()  
mybrusfig.savefig('signal_strength_vs_noise.jpg')


fig_all.suptitle('Mean signal as function of temperature')
fig_all.supxlabel('Temperature [C]')
fig_all.savefig('FlatfieldTempDependence.jpg')
fig_all.legend()


