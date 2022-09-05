#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 09:33:32 2020

@author: lindamegner
"""


    

from mats_l1_processing.L1_calibration_functions import compensate_bad_columns #,  get_true_image
from mats_l1_processing.read_in_functions import read_all_files_in_directory
import matplotlib.pyplot as plt
import numpy as np
from mats_l1_processing.LindasCalibrationFunctions import  plot_CCDimage   


def columnbin(image,binsize):
    orgsize=image.shape
    newsize=[orgsize[0],int(orgsize[1]/binsize)]
    binnedimage=np.ones(newsize)
    for i in range(0,newsize[1]):
        binnedimage[:,i]=image[:,i*binsize:i*binsize+binsize].sum(1)
    return binnedimage
 
def simbin(image,tblnk,binsize):
    blank_off=tblnk-128
    image=image-blank_off-128
    simbinimage=columnbin(image,binsize)+blank_off+128
    return simbinimage


def diffplot2(image1, image2, title1, title2, fig, ax,axind, clim=999, climdiff=999):
    from mats_l1_processing.LindasCalibrationFunctions import  plot_CCDimage    
    diffimg=image1-image2
    if clim==999:        
        plot_CCDimage(image1,fig, ax[0,axind], title=title1)
        plot_CCDimage(image2,fig, ax[1,axind], title=title2)
    else:
        plot_CCDimage(image1,fig, ax[0,axind], title=title1, clim=clim)
        plot_CCDimage(image2,fig, ax[1,axind], title=title2, clim=clim)

    if climdiff==999:   
        plot_CCDimage(diffimg,fig, ax[2,axind], title='row 1 - row 2')
    else:
        plot_CCDimage(diffimg,fig, ax[2,axind], title='row 1 - row 2', clim=climdiff)
        
    ax[0,axind].set_aspect('auto')
    ax[1,axind].set_aspect('auto')
    ax[2,axind].set_aspect('auto')
    return fig



def doublediffplot(fig,ax,image1, image2, title1, title2, clim=999, climdiff=999, icol=0):
    from mats_l1_processing.LindasCalibrationFunctions import  plot_CCDimage    

    diffimg=image1-image2
    if clim==999:        
        plot_CCDimage(image1,fig, ax[0,icol], title=title1)
        plot_CCDimage(image2,fig, ax[1,icol], title=title2)
    else:
        plot_CCDimage(image1,fig, ax[0,icol], title=title1, clim=clim)
        plot_CCDimage(image2,fig, ax[1,icol], title=title2, clim=clim)
    if climdiff==999:   
        plot_CCDimage(diffimg,fig, ax[2,icol], title='pic 1-pic2')
    else:
        plot_CCDimage(diffimg,fig, ax[2,icol], title='pic 1-pic2', clim=climdiff)
        

    return

def read_all_files_in_RacFiles_out(directory):
    CCDitems=read_all_files_in_directory('rac',directory+'RacFiles_out/')
    return CCDitems
    

directories=[]
directories.append('/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/FunctionalTests/bad_column_test/sync_2/' )
directories.append('/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/FunctionalTests/bad_column_test/sync_0/' )


figax=[]
fig=[]
ax=[]
fig1, ax1 =plt.subplots(3,2)
image_binned_no_bc_sync2=[]
comp_image_L1_sync2=[]
for idir, directory in enumerate(directories):
    CCDitems=read_all_files_in_RacFiles_out(directory)
    
    
    item_nonbinned=CCDitems[0]
    image_nonbinned=CCDitems[0]['IMAGE']
    binsize=10
    image_selfbinned=simbin(image_nonbinned, item_nonbinned['TBLNK'],binsize)
    image_binned_no_bc=CCDitems[1]['IMAGE']
    
    plotfirstfigure=True
    if plotfirstfigure: 
        clim=[1070,1280]

        doublediffplot(fig1,ax1,image_binned_no_bc, image_selfbinned, title1='Binned on chip', title2='Binned in postprocessing', clim=clim, icol=idir)
        #fig.savefig('Bintest.jpg')
        for iax in fig1.get_axes():
            iax.set_aspect('auto')

        #fig1.get_axes()[5].set_aspect('auto')
    
    clim=[800,1300]
    itestexample=5

    for i, CCDitem in enumerate(CCDitems[2:]):  

        print('Bad columns:',CCDitem['BC'] )
        #fig=diffplot2(image_binned, CCDitem['IMAGE'],'binned no bc','non-compensated', fig, ax,axind=0, clim=clim)

        if directory=='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/FunctionalTests/bad_column_test/sync_2/':
            sync='SYNC: 2'
            comp_image=compensate_bad_columns(CCDitem)
            figax.append(plt.subplots(3,4, figsize=(14,12)))
            fig.append(figax[i][0])
            ax.append(figax[i][1])
            # Plot the non-compensated
            fig[i]=diffplot2(image_binned_no_bc, CCDitem['IMAGE'],'binned no bc '+sync,'bc uncompensated '+sync,fig[i], ax[i],axind=0, clim=clim)      
            # Plot the L1 Georgi algorithm compensated 
            fig[i]=diffplot2(image_binned_no_bc, comp_image,'binned no bc '+sync,'post-compensated '+sync, fig[i], ax[i],axind=1, clim=clim)
            image_binned_no_bc_sync2.append(image_binned_no_bc)
            comp_image_L1_sync2.append(comp_image)
            
            if i==itestexample: 
                image_uncompensated=CCDitem['IMAGE']
        else:
            sync='SYNC: 0'         
            # Plot the OBC-compensated
            fig[i]=diffplot2(image_binned_no_bc, CCDitem['IMAGE'],'binned no bc '+sync,'OBC compensated '+sync, fig[i], ax[i],axind=2, clim=clim)
 
             # Plot the 4th column, the difference between column 2 and 3
            fig[i]=diffplot2(image_binned_no_bc_sync2[i]-image_binned_no_bc,comp_image_L1_sync2[i]-CCDitem['IMAGE'] ,'col 2 - col 3 ','col 2 - col 3', fig[i], ax[i],axind=3, clim=[-100,100])
            
            
            if i==itestexample: 
                image_OBCcompensated=CCDitem['IMAGE']
            # txtstring1='mean = '+str(np.mean(image_binned_no_bc_sync2[i]-image_binned_no_bc))
            # txtstring2='std = '+str(round(np.std(image_binned_no_bc_sync2[i]-image_binned_no_bc),2))
            # ax[i][0,3].text(0.5, 10,txtstring1)
            # ax[i][0,3].text(0.5, 20,txtstring2)
            
        fig[i].suptitle('Test: '+ str(i)+' Bad columns:'+str(CCDitem['BC']))
        
        fig[i].savefig('BadColumn_Compensated_OBC_vs_after'+str(i)+'.jpg')
        

            

        
    fig1.savefig('Binned_on_chip_versus_using_after_processing.jpg')
    

# Make example plot for paper:
paperfig, paperax =plt.subplots(2,1)
plot_CCDimage(image_uncompensated,paperfig, paperax[0], title='', clim=[900, 1300])
plot_CCDimage(image_OBCcompensated,paperfig, paperax[1], title='', clim=[900, 1300])
paperfig.savefig('Uncompensated_vs_compensated.jpg')

#image_bias_sub = get_true_image(CCDitem)