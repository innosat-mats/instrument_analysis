#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 12:49:15 2021

@author: lindamegner
"""


from mats_l1_processing.items_units_functions import read_files_in_protocol_as_ItemsUnits
from mats_l1_processing.experimental_utils import readprotocol
from matplotlib import pyplot as plt
import numpy as np



def findpattern(lena,template):
    from scipy import signal
    orgimage=lena.copy()
    lena = lena-lena.mean()
    template= template- template.mean()
    corr = signal.correlate2d(lena, template, boundary='fill', mode='full')
    y, x = np.unravel_index(np.argmax(corr), corr.shape) # find the match    
    bestcorrpic=orgimage[y-template.shape[0]+1:y+1,x-template.shape[1]+1:x+1]




    # res_h_pic=bestcorrpic[54:62,50:80]
    # res_v_pic=bestcorrpic[45:70,85:93]
    # res_h=np.mean(res_h_pic,0)
    # res_v=np.mean(res_v_pic,1)
    # bright_pic=bestcorrpic[20:40,5:20]
    # bright=np.mean(bright_pic)
    # dark_pic=bestcorrpic[0:10,0:10]  
    # dark=np.mean(dark_pic)  
    

    return x , y, bestcorrpic







output_imagedir='images/'

directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210215OHBLimbImage/'
df=readprotocol(directory+'protocol_dark_bright_100um_incl_IR3.txt')

read_from='rac'





fig, axes = plt.subplots(6,3,figsize=(10,10))

ItemsUnits=read_files_in_protocol_as_ItemsUnits(df,directory,3, read_from)


ind=0
template=ItemsUnits[ind].subpic[230:270,850:890]



#tempfig, tempaxes = plt.subplots(1)
#sp1=plot_CCDimage(template,tempfig, tempaxes,title='template'+ItemsUnits[ind].imageItem['channel'])

#fig1, axes1 = plt.subplots(6,1,figsize=(10,3))
fig1 = plt.figure(figsize=(5,7))

for i, ItemsUnit in enumerate(ItemsUnits[:]):



    for j in range(0,3):
        if j==2:  clim=[0, 50]
        else: clim=999
        ItemsUnit.plot(fig,axes[i,j], whichpic=j, title=ItemsUnit.imageItem['channel'], clim=clim)    
    if ItemsUnit.imageItem['channel']=='IR2' or ItemsUnit.imageItem['channel']=='IR4':
        image=np.fliplr(ItemsUnit.subpic)
    else:    
        image=ItemsUnit.subpic
    
    x, y, bestcorrpic= findpattern(image,template)


#    figLarge, axesLarge = plt.subplots(1)    
   # plot_CCDimage(image, figLarge, axesLarge, title='Image '+ItemsUnit.imageItem['channel'])
    
    
    ycenter, xcenter=np.unravel_index(bestcorrpic.argmax(),bestcorrpic.shape)

    
    if ItemsUnit.imageItem['channel']=='IR1' or ItemsUnit.imageItem['channel']=='IR2' :
        vres=0.4/2 #km
        hres=10/2 # km
    elif ItemsUnit.imageItem['channel']=='IR3' or ItemsUnit.imageItem['channel']=='IR4' :
        vres=0.8/2 #km
        hres=50/2 # km       
    elif ItemsUnit.imageItem['channel']=='UV1'or ItemsUnit.imageItem['channel']=='UV2' :
        vres=0.2/2 #km
        hres= 5/2 #km





    hres=hres/10
    xpix=hres/0.13 #konvert to pixels
    ypix=vres/0.13
    
    mask = np.zeros(bestcorrpic.shape, dtype=bool)
    ymin=int(ycenter-ypix)
    ymax=int(ycenter+ypix)
    xmin=int(xcenter-xpix)
    xmax=int(xcenter+xpix)
    

    if ymin<0: ymin=0
    if ymax>mask.shape[0]: ymax=mask.shape[0]
    if xmin<0: xmin=0
    if xmax>mask.shape[1]: xmax=mask.shape[1]    
    
    # 
    # mask[ymin:ymax+1,xmin:xmax+1] = True
    # largemask = np.zeros(bestcorrpic.shape, dtype=bool)
    # largemask[int(ycenter-4*ypix):int(ycenter+4*ypix),int(xcenter-2*xpix):int(xcenter+2*xpix)] = True
    
    # Exclude x: 
    mask[ymin:ymax+1,:] = True

    largemask = np.zeros(bestcorrpic.shape, dtype=bool)
    largemask[int(ycenter-4*ypix):int(ycenter+4*ypix),:] = True
    
# # Circular or ellipsoid masks below:
#     ellipse=False
#     if ellipse:
#         xrad=hres/0.13 #konvert to pixels
#         yrad=vres/0.13
#     else:#circle
#         xrad=vres/0.13
#         yrad=vres/0.13
#     angle = np.linspace( 0 , 2 * np.pi , 150 ) 
#     xval = xrad * np.cos( angle ) 
#     yval = yrad * np.sin( angle ) 
    
#     xlen, ylen=bestcorrpic.shape
#     xmask = np.arange(0, xlen)
#     ymask= np.arange(0, ylen)    
#     mask = (xmask[np.newaxis,:]-xcenter)**2/xrad**2 + (ymask[:,np.newaxis]-ycenter)**2/yrad**2 < 1**2
#     largemask = (xmask[np.newaxis,:]-xcenter)**2/xrad**2 + (ymask[:,np.newaxis]-ycenter)**2/yrad**2 < 4**2
    
    #Calibrate to values outside the large mask
    bestcorrpic=bestcorrpic-np.average(bestcorrpic[~largemask]) #remove background
    insidesum=np.sum(bestcorrpic[mask])
    large_circle_sum=np.sum(bestcorrpic[largemask])
    wholesum=np.sum(bestcorrpic)
    fraction=round(insidesum/large_circle_sum,2)
    #bestcorrpic[mask]=9999


    #axes1[i,0].imshow(image[150:350, 800:1200], cmap=plt.cm.gray)   
    #axes1[i,0].set_title('Image '+ItemsUnit.imageItem['channel'])

    axes1 = fig1.add_subplot(3, 2, i+1)
    axes1.pcolormesh(bestcorrpic)   
    axes1.plot([0,bestcorrpic.shape[1]-1], [ymin,ymin],'r')
    axes1.plot([0,bestcorrpic.shape[1]-1], [ymax,ymax],'r')
    axes1.set_title(ItemsUnit.imageItem['channel']+' Part of signal:'+str(fraction))
    #axes1.plot( xcenter+xval, ycenter+yval, color='r' ) 
    #axes1[i].text(2,2,'part of signal:'+str(insidesum/large_circle_sum), color='r')
        
    #print('sum: ',insidesum, 'part of tot signal: ', insidesum/wholesum )    
    
    
    fig1.savefig(output_imagedir+'LimbImage100um.jpg')
    
    
    
    
