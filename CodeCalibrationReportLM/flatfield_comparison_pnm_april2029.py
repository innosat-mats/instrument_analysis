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

#from LindasCalibrationFunctions import plot_CCDimage


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

channellist=['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']
channellist=['IR3', 'UV2']
for channel in channellist:

        
    # #dir1='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/PayloadImages/'
    # dir1='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/PayloadImages20210421_16_18/'
    # if channel=='IR1':
    #     filename1=dir1+'1303052550753753600_1'+'.pnm'
    #     dfilename1=dir1+'1303052736716888320_1'+'.pnm'
    # elif channel=='IR2':
    #     filename1=dir1+'1303053042442397952_4'+'.pnm'
    #     dfilename1=dir1+'1303053108130874624_4'+'.pnm'
    # elif channel=='IR3':
    #     filename1=dir1+'1303053614869430528_3'+'.pnm'
    #     dfilename1=dir1+'1303053680075119104_3'+'.pnm'
    # elif channel=='IR4':
    #     filename1=dir1+'1303053921036712704_2'+'.pnm'
    #     dfilename1=dir1+'1303053986963058432_2'+'.pnm'
    # elif channel=='UV1':
    #     filename1=dir1+'1303054301145446656_6'+'.pnm'
    #     dfilename1=dir1+'1303054436714050304_6'+'.pnm'
    # elif channel=='UV2':
    #     filename1=dir1+'1303055407319107072_5'+'.pnm'
    #     dfilename1=dir1+'1303055745598846464_5'+'.pnm'    
    
   
    dir1='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_Omnisys/ScreenCalibration/'
    if channel=='IR3':
        #filename1=dir1+'rawoutput'+'12021'+'.pnm' 
        #dfilename1=dir1+'rawoutput'+'14198'+'.pnm'   
        filename1=dir1+'output'+'3778'+'.pnm' 
        dfilename1=dir1+'output'+'757'+'.pnm'  
    elif channel=='UV2':
        filename1=dir1+'rawoutput'+'985'+'.pnm' #Jörg
        dfilename1=dir1+'rawoutput'+'3243'+'.pnm'     
    else:
        print('Error no',channel,' channel')
     
    
    
    when='May2020'
    
    if when=='May2020':
        dir2='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200506_flatfields_roomtemp/PayloadImages/'
        if channel=='IR1':
            filename2=dir2+'3713_10412_1'+'.pnm' #Jörg
            dfilename2=dir2+'3790_54246_1'+'.pnm' #Jörg
        elif channel=='IR2':
            filename2=dir2+'10501_10389_4'+'.pnm'
            dfilename2=dir2+'10574_55978_4'+'.pnm'
        elif channel=='IR3':
            filename2=dir2+'12414_4234_3'+'.pnm'
            dfilename2=dir2+'12485_41030_3'+'.pnm'
        elif channel=='IR4':
            filename2=dir2+'14775_11476_2'+'.pnm'
            dfilename2=dir2+'14836_31942_2'+'.pnm'
        elif channel=='UV1':
            filename2=dir2+'22264_47961_5'+'.pnm'
            dfilename2=dir2+'22463_26703_5'+'.pnm'
        elif channel=='UV2':
            filename2=dir2+'18430_43035_6'+'.pnm'
            dfilename2=dir2+'18638_6404_6'+'.pnm'
            
    elif when=='May2019':        
        dir2='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_AIT/20190520/PayloadImages/'
        if channel=='IR1':
            filename2=dir2+'rawoutput'+'11220'+'.pnm' #Jörg
            dfilename2=dir2+'rawoutput'+'9126'+'.pnm' 
        elif channel=='IR2':
            filename2=dir2+'rawoutput'+'15462'+'.pnm'
            dfilename2=dir2+'rawoutput'+'13362'+'.pnm'
        elif channel=='IR3':
            filename2=dir2+'rawoutput'+'1278'+'.pnm'
            dfilename2=dir2+'rawoutput'+'636'+'.pnm'
        elif channel=='IR4':
            filename2=dir2+'rawoutput'+'5503'+'.pnm'
            dfilename2=dir2+'rawoutput'+'7002'+'.pnm'
        elif channel=='UV1':
            filename2=dir2+'rawoutput'+'7368'+'.pnm'
            dfilename2=dir2+'rawoutput'+'2313'+'.pnm'
        elif channel=='UV2':
            filename2=dir2+'rawoutput'+'4043'+'.pnm'
            dfilename2=dir2+'rawoutput'+'7803'+'.pnm'
            
     

    print('file', filename2)
    fig,ax=plt.subplots(3,1,figsize=(6,5.5))
    pic1 = np.float64(Image.open(filename1))  # read image
    pic2 = np.float64(Image.open(filename2))
    
    picd1 = np.float64(Image.open(dfilename1))  # read image
    picd2 = np.float64(Image.open(dfilename2))
    
    img1=pic1-picd1
    img2=(pic2-picd2)
    img2_scaled=img2*img1.mean()/img2.mean()
    diff=img1-img2_scaled
    
    
    plot_CCDimage(img1, fig, ax[0],title='Omnisys April 2019 B-D')
    plot_CCDimage(img2_scaled, fig, ax[1],title=when+', B-D scaled to same mean')
    plot_CCDimage(diff, fig, ax[2],title='difference panel1-panel2')
    
    fig.suptitle(channel)
    fig.savefig('Limb_flatfield_diff_Omnisys2019_to_'+when+ channel+'.jpg')