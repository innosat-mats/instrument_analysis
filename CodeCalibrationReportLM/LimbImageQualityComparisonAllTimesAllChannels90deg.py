#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:00:16 2021

@author: lindamegner

Routine to compare image qualtiy at different times. 

This routine does not subract the dark current picutre since there are no dark current pictures in the older data.

"""


from mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage, read_all_files_in_protocol 
from mats_l1_processing.read_in_functions import readprotocol
import matplotlib.pyplot as plt
import numpy as np
#from jpeglib import read12bit_jpegfile
from PIL import Image

def findpattern(lena,template):
    from scipy import signal
    orgimage=lena.copy()
    lena = lena-lena.mean()
    template= template- template.mean()
    corr = signal.correlate2d(lena, template, boundary='fill', mode='full')
    y, x = np.unravel_index(np.argmax(corr), corr.shape) # find the match    
    bestcorrpic=orgimage[y-template.shape[0]+1:y+1,x-template.shape[1]+1:x+1]




    res_h_pic=bestcorrpic[54:62,50:80]
    res_v_pic=bestcorrpic[45:70,85:93]
    res_h=np.mean(res_h_pic,0)
    res_v=np.mean(res_v_pic,1)
    bright_pic=bestcorrpic[20:40,5:20]
    bright=np.mean(bright_pic)
    dark_pic=bestcorrpic[0:4,0:20]  
    dark=np.mean(dark_pic)  
    
    
    
    if True:
        fig, axes = plt.subplots(3,2)
        s=axes[2,0].imshow(res_h_pic)   
        axes[2,0].set_title('res_h_pic')
        fig.colorbar(s, ax=axes[2,0])
        s=axes[2,1].imshow(res_v_pic)   
        axes[2,1].set_title('res_v_pic')
        fig.colorbar(s, ax=axes[2,1])
        s=axes[1,0].imshow(bright_pic)   
        axes[1,0].set_title('bright_pic')
        fig.colorbar(s, ax=axes[1,0])
        s=axes[1,1].imshow(dark_pic)   
        axes[1,1].set_title('dark_pic')   
        fig.colorbar(s, ax=axes[1,1])
        s=axes[0,0].imshow(lena)   
        axes[0,0].set_title('original picture')
        fig.colorbar(s, ax=axes[0,0])
        s=axes[0,1].imshow(bestcorrpic)   
        axes[0,1].set_title('pattern selected area')
        fig.colorbar(s, ax=axes[0,1])


        

    
    return x , y, bestcorrpic, res_h, res_v, bright, dark

output_imagedir='images/'
#Setup directories 
directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210215OHBLimbImage/'
df_protocol=readprotocol(directory+'protocol_dark_bright_90degree_incl_IR3.txt')
read_from='rac'  
df_bright=df_protocol[df_protocol.DarkBright=='B']
CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)


woojindir='/Users/lindamegner/MATS/retrieval/Calibration/USAF_raw_Woojin_190321-190326/'
imagedir='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_AIT/20190517/data/'
arviddir='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_Omnisys/20190425Arvid/'
afterglue='/Users/lindamegner/MATS/retrieval/Calibration/FM_tests_after_glue/20190918/'
cooltestdir22='/Users/lindamegner/MATS/retrieval/Calibration/CoolingTests/PayloadImages20200114/'
cooltestdir14='/Users/lindamegner/MATS/retrieval/Calibration/CoolingTests/PayloadImages20200109/'
cooltestdir11='/Users/lindamegner/MATS/retrieval/Calibration/CoolingTests/PayloadImages20200109/'
cooltestdir05='/Users/lindamegner/MATS/retrieval/Calibration/CoolingTests/PayloadImages20200110/'


# =============================================================================
# Names on Box:
#     PayloadImages20200114 is called 20200114
#     PayloadImages20200109 is called 20200109
#     PayloadImages20200110 is called 20200110
# =============================================================================





fig3, axes3 = plt.subplots(5,2,figsize=(10,10))
channellist=['IR1', 'IR2','IR3', 'IR4','UV2']
#channellist=['IR3']

#Select channel and choos picture for that channel

for ichannel, channel in enumerate(channellist):


    
    CCDitem_select= list(filter(lambda x: ( x['channel']==channel),CCDitems))
    CCDitem=CCDitem_select[0]
    
    
    
    
    
    if channel == 'IR1':
        arvidpicID=[6590,6725]
        woojinpicID=14260
        pics=[3179, 5009, 4082]
        gpics=[3604]
        #cool15=['2328_58444_1','2527_49307_1'] 
        cool22=['1112_62821_1','1258_44467_1.']      
        cool14=['76976_5626_1', '77165_18071_1']  
        cool11=['86746_37588_1']
        cool05=['165752_3478_1','165906_59738_1 ']
    elif channel == 'IR2':
        arvidpicID=[6342,5610]
        woojinpicID=2595
        pics=[5513, 5677, 6329]
        gpics=[4485]
        #cool15=['2926_47785_4','2817_47370_4'] 
        cool22=['1608_10394_4','1718_4608_4']      
        cool14=['77454_20269_4', '77550_18550_4']
        cool11=['87070_48650_4']
        cool05=['166344_9175_4','166417_2165_4']
    elif channel == 'IR3':
        arvidpicID=[10622,10504]
        woojinpicID=479
        pics=[7516,7644, 7974]
        gpics=[4707]
        #cool15=['1595_54355_3','1844_53355_3'] 
        cool22=['1984_57720_3','2077_65069_3']      
        cool14=['77842_7626_3','77964_18679_3']
        cool11=['87511_47392_3']
        cool05=['166611_4717_3','166677_6728_3']
    elif channel == 'IR4':
        arvidpicID=[10912,11048]
        woojinpicID=8584
        pics=[8199,8440,8619]
        gpics=[5054]    
        cool15=['2055_40291_2','2146_52709_2']
        cool22=['2307_1620_2','2390_40237_2']  
        cool11=['87841_42996_2']
        cool05=['167103_11332_2','167181_5036_2']  
        cool14=['78224_19835_2','78330_18412_2']
    elif channel == 'UV2':
        arvidpicID=[12685,11958]
        woojinpicID=13307
        pics=[9164,9643,10252]
        gpics=[5730]
        #cool15=['3538_51913_6','3293_48923_6']
        cool22=['2624_58115_6','2700_41090_6']    
        cool14=['78622_18653_6','78702_7393_6']
        cool05=['167482_44_6','167686_870_6']
        cool11=['88153_51667_6']
    elif channel == 'UV1':
        #cool15=['','']
        cool22=['3318_11356_5','3636_37757_5']
        cool22d=['2964_40791_5','3953_44638_5']
        cool14=['79512_10236_5','79845_13180_5']
        cool14d=['79179_23973_5','80252_21365_5']
        cool11=['88679_59108_5','88976_19080_5']
        cool11d=['88370_48899_5','89276_62313_5']
        cool05=['168076_32759_5','168388_62696_5']
        cool05d=['167777_42868_5','168708_51846_5']
    
    
    
    
    
    ############
    target=0 # 0= USAF 1=FOV
    targetname=['USAF', '100um']
    
    # woojinpic=2**4*(np.int64(read12bit_jpegfile(woojindir +channel + '/'+'USAF__output'+str(woojinpicID)+'.jpg')))
    # arvidpic=2**4*np.int64(read12bit_jpegfile(arviddir + channel + '/'+'output'+str(arvidpicID[target])+'.jpg'))
    # AITpic=2**4*np.flip(np.int64(read12bit_jpegfile(imagedir + 'output'+str(pics[2])+'.jpg')))
    # AITpic_after_glue=2**4*np.flip(np.int64(read12bit_jpegfile(afterglue + 'output'+str(gpics[0])+'.jpg')))
    
    woojinpic=2**4*(np.int64(Image.open(woojindir +channel + '/'+'USAF__output'+str(woojinpicID)+'.pnm')))
    arvidpic=2**4*np.int64(Image.open(arviddir + channel + '/'+'USAF__output'+str(arvidpicID[target])+'.pnm'))
    AITpic=2**4*np.flip(np.int64(Image.open(imagedir + 'output'+str(pics[2])+'.pnm')))
    AITpic_after_glue=2**4*np.flip(np.int64(Image.open(afterglue + 'output'+str(gpics[0])+'.pnm')))
    
    
    
    
    cooltest_22C=np.int64(Image.open(cooltestdir22+cool22[target]+'.pnm'))
    cooltest_14C=np.int64(Image.open(cooltestdir14+cool14[target]+'.pnm'))
    cooltest_11C=np.int64(Image.open(cooltestdir11+cool11[target]+'.pnm'))
    cooltest_05C=np.int64(Image.open(cooltestdir05+cool05[target]+'.pnm'))
    
    
    
    
    
    # Use woojin IR2 as template because that is very sharp
    
    templatepic=2**4*(np.int64(Image.open(woojindir +'IR2' + '/'+'USAF__output'+str(2595)+'.pnm')))
    template=np.fliplr(np.copy(templatepic[170:244,820:940]))
    #template=np.rot90(template)
    
    
    
    tempfig, tempaxes = plt.subplots(1)
    sp1=plot_CCDimage(template,tempfig, tempaxes,title='template Woojin IR2')  
    
    
    
    
    
    
    piclist=[woojinpic, arvidpic, AITpic, AITpic_after_glue, cooltest_22C, cooltest_14C, cooltest_11C, cooltest_05C, CCDitem['IMAGE']]
    namelist=['Marsh 2018','April 2018','May 2019 ', 'Sept 2019', 'Jan 2020 22C','Jan 2020 14C','Jan 2020 11C','Jan 2020 5 C', 'Feb 2021' ]
    rotated=[False, True, True, True, True, True, True, True, True]
    picdictlist=[]
    for index, pic in  enumerate(piclist):
        picdict= {
            "pic": pic,
            "name": namelist[index],
            "rotated": rotated[index]
            }
        picdictlist.append(picdict)    
    picdictlistselect=[picdictlist[i] for i in [4,8]]
    #picdictlistselect=picdictlist
    
    fig, axes = plt.subplots(len(picdictlistselect),1,figsize=(14,10))
    
    for index, picdict in enumerate(picdictlistselect):
        plot_CCDimage(picdict['pic'],fig, axes[index],title=namelist[index]+' '+channel)
    
    
    
    resvec_h=[]
    resvec_v=[]
    fig2, axes2 = plt.subplots(len(picdictlistselect),1,figsize=(14,10))
    #fig3, axes3 = plt.subplots(2,1,figsize=(10,10))
        
    
    for index, picdict in enumerate(picdictlistselect):
        pic=picdict['pic']
        if channel=='IR2' or channel=='IR4':
            pic=np.fliplr(pic)
            
    
            
        if picdict['rotated']:    
            pic_input=np.rot90(np.copy(pic), 3)    # rotates 3*90 degrees
        else:
            pic_input=np.copy(pic)
        x, y, bestcorrpic, res_h, res_v, bright, dark= findpattern(pic_input,template)
        
        if picdict['rotated']: #swich horizontla and vertical resolution
            res_v_temp=res_v
            res_v=res_h
            res_h=res_v_temp
            
     
        sp=axes3[ichannel, 0].plot((res_h-dark)/(bright-dark),label=picdict['name'])
        axes3[ichannel, 0].set_title('Horizonal resolution '+channel)
        axes3[ichannel, 0].legend()
        sp=axes3[ichannel, 1].plot((res_v-dark)/(bright-dark),label=picdict['name'])   
        axes3[ichannel, 1].set_title('Vertical resolution '+channel)
        axes3[ichannel, 1].legend()
        if picdict['rotated']:
            plot_CCDimage(bestcorrpic,fig2,axes2[index], title=channel+' '+picdict['name']+' rotated 90 degree')
        else:
            plot_CCDimage(bestcorrpic,fig2,axes2[index], title=channel+' '+picdict['name'])
            
        resvec_h.append(res_h)
        resvec_v.append(res_v)

    #fig2.savefig('Limb_image_quality_alltime'+channel+'.jpg')    
fig3.savefig(output_imagedir+'image_resolution_alltime_all_channels_90deg.jpg')
