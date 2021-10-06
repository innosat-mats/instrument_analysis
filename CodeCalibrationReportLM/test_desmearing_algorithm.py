#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:19:08 2020

@author: lindamegner


This script is currently not running /LM20210209
"""


from mats_l1_processing.read_in_functions import read_CCDitems
from database_generation.read_in_imgview_functions import readselectedimageviewpics


from mats_l1_processing.get_temperature import create_temperature_info_array, add_temperature_info

from mats_l1_processing.L1_calibration_functions import desmear_true_image
import datetime 
from cali_ReadJorgsFiles import readjorgfile
from mats_l1_processing.LindasCalibrationFunctions import UniqueValuesInKey
from mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage, plot_CCDimage_hmean

import numpy as np
import matplotlib.pyplot as plt



# def plot(image,title,darkarea):
#     plt.figure()
#     mean_img=image.mean()
#     mean_at_dark=np.mean(image[darkarea[0]:darkarea[1],darkarea[2]:darkarea[3]])
# #    print('mean at dark=', mean_at_dark)
#     std=image.std()
#     p0= plt.imshow(image, aspect='auto',origin='lower', vmin=mean_img-2*std, vmax=mean_img+2*std)
#     plt.title(title+ ' mean= '+str(mean_img))
#     plt.xlabel('Pixels')    
#     plt.ylabel('Pixels')
#     plt.colorbar()
#     plt.text(150,100,'mean at dark= '+ str(mean_at_dark))

#     return  p0

# def CCDplot(fig,axis,pic,title='',clim=999):
                    
#     sp=axis.pcolormesh(pic)
#     axis.set_title(title)
#     if clim==999:
#         mean=pic.mean()
#         std=pic.std()
#         sp.set_clim([mean-2*std,mean+2*std])
#     else:
#         sp.set_clim(clim)

#     fig.colorbar(sp,ax=axis)

#     return sp

def CCDplot(fig,axis,image,title='',clim=999):
    sp=plot_CCDimage(image, fig, axis, title=title)
    return sp


    
    
def meanandstd(image):
    mean=np.mean(image)
    std=np.std(image)
    return mean, std


        
        

#     #################################################
#     # Read in the data                            #
#     #################################################
    
#      # 0 read_from=0 is from KTH images, 1 is fron rac file, 2 is from image viewer
    
#     if read_from==0: #Read KTH files
#         CCDitem, flag =readimage_create_CCDitem('/Users/lindamegner/MATS/retrieval/Level1/data/2019-02-08 rand6/', 1)
#     elif read_from==1: #Read from rac file new version as of November 2019 
#         rac_image_json_file='rac20191106testme/images.json'
#         rac_packets_json_file='rac20191106testme/packets.json'    
#     #    rac_image_json_file='rac20190818-152721/images.json'
#         rac_image_json_dir='/Users/lindamegner/MATS/retrieval/Level0/MATS-L0-processing-master/'
#         CCDitem=read_CCDitems(rac_image_json_file,rac_image_json_dir)
#     elif read_from==2: #read image and textfile created by image viewer
#         rawflag=1
#         dirname='/Users/lindamegner/MATS/retrieval/Calibration/FM_tests_after_glue/20191106/'
#         picnr=14976
#         CCDitem=readimageviewpic(dirname,picnr,rawflag)
        
#     return CCDitem

SigM=1   

channel='UV2'
axes=[None]*2
figs=[None]*2
figind=0


NrOfExpTimes=7
figs[0], axes[0] = plt.subplots(NrOfExpTimes,4,figsize=(10,8))
figs[0].set_size_inches(10,10)
figs[0].suptitle(str(channel)+' SigM='+str(SigM)+ ', B-D, scaled to 1s exptime ')


figs[1], axes[1] = plt.subplots(NrOfExpTimes,4,figsize=(10,8))
figs[1].set_size_inches(10,10)
figs[1].suptitle(str(channel)+' SigM='+str(SigM)+ ', B-D for no shutter desmeared (left) with shutter(mid), difference (right), scaled to 1s exp time')


CCD_noshutter=[]
CCD_shutter=[]

read_from='imgview'
epoch=datetime.datetime(1980,1,6)
labtemp=3 #Currently not used.
shuttertime=2 #time to compensate the exposure time for that the shutter is open less time (normally 2s)




picdirnames=['/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200322_flatfields_MISU/',
             '/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200330_flatfields_MISU/']
# =============================================================================
# NOTE!!! The second directory above is the same as 
# '/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200330_flatfields_0C/PayloadImagesOld']
#
#  Names on box: 
#  
#  The first directory is called 20200322_flatfields_MISU_no_shutter
#  The second directory is called 20200330_flatfields_MISU on box which is the same as directory 
#  20200330_flatfields_0C
#
# =============================================================================

    

if SigM==1:
    jorgtxtfiles=['flatfields_200322_SigMod1.txt','flatfields_200330_SigMod1.txt']
elif SigM==0:
    jorgtxtfiles=['flatfields_200322_SigMod0.txt', 'flatfields_200330_SigMod0.txt'] 

for ifile, picdirname in enumerate(picdirnames): 
    jorgtxtfile=jorgtxtfiles[ifile]
    df=readjorgfile(picdirname+'/'+jorgtxtfile)
    
    
    rac_image_json_file='images.json'
    rac_packets_json_file='packets.json'     
    #    rac_image_json_file='rac20190818-152721/images.json'
    rac_sub_dir='rac20191106testme/'
    retrieval_dir='/Users/lindamegner/MATS/retrieval/Level0/MATS-L0-processing-master/'
    
    
    
      
    if read_from=='rac':
        CCDitems=read_CCDitems(rac_image_json_file,rac_sub_dir, retrieval_dir)
    elif read_from=='imgview': #read image and textfile created by image viewer
    #    rawflag=1
    
       # CCDitems=readimageviewpics(dirname,rawflag)
        CCDitems=readselectedimageviewpics(picdirname,list(df['ID']))
    
    
    # Create temperature information array. 
    #This should be done once for every json packet
    if read_from=='rac':    
        temperaturedata, relativetimedata=create_temperature_info_array(retrieval_dir+rac_sub_dir+rac_packets_json_file)
    else:
        temperaturedata=999
        relativetimedata=999
    
    for CCDitem in CCDitems:
        CCDitem['reltime']=int(CCDitem['EXPTS'])+int(CCDitem['EXPTSS'])/2**16 
        CCDitem['read_from']=read_from
        CCDitem=add_temperature_info(CCDitem,temperaturedata,relativetimedata,labtemp)
        timestamp=epoch+datetime.timedelta(0,CCDitem['reltime'])
        
    #    print(timestamp)
    #    print(CCDitem['temperature'])
        
    
    
    
    #Add dict entry that tells if it is a dark picture or an image
    for CCDitem in CCDitems:   
    #tmpind2
    #    tmpindex=df.loc[df['ID']==CCDitem['ID']].index[0]

        tmpind=df.index[df.ID == CCDitem['id']]
    
        if tmpind.size==0 : 
             CCDitem['ImgOrDark']=-1
        elif tmpind.size==1:
             CCDitem['ImgOrDark']=int(df['ImgOrDark'][tmpind[0]])
        elif tmpind.size>1:
             CCDitem['ImgOrDark']=999   
             print('Warning: Two of the same picture ID in Jörgs table')
    
    
    #Sort list and devide in channels:
             
    CCDitems.sort(key=lambda x: x['channel'])
        
    
      
    
    #Now select CCD pictures to analyse   
    
    CCDitems_channel = list(filter(lambda x: (x['channel']==channel),CCDitems))
    # if channel=='IR1':
    #      imclim=700 
    #      sublim=60
    # elif channel=='IR2':
    #      imclim=2000
    #      sublim=100     
    # elif channel=='IR3':
    #      imclim=300
    #      sublim=60     
    # elif channel=='IR4':
    #      imclim=300
    #      sublim=60     
    # elif channel=='UV1':
    #      imclim=30
    #      sublim=5     
    # elif channel=='UV2':
    #      imclim=100
    #      sublim=40     
    
    #CCDarray=np.array(CCDitems_channel)
    
    
    
    
    #CCDitems_channel.sort(key=lambda x: x['TEXPMS'])
    
    
    ExpTimes=UniqueValuesInKey(CCDitems_channel,'TEXPMS')
    ExpTimes[:] = [int(x)/1000 for x in ExpTimes if x/1000>shuttertime]
    
    
    
    print('ifile ', jorgtxtfile, ' Exposure times from rac ', ExpTimes)
    
    
    
    list_scaled=[]
    if ifile==0:
        ExpTimes=[x for x in ExpTimes[:NrOfExpTimes]]
        ExpTimesFirstFile=ExpTimes
    elif ifile==1:
        start=ExpTimes.index(ExpTimesFirstFile[0]+shuttertime)
        ExpTimes=[x for x in ExpTimes[start:start+NrOfExpTimes]]
        if start ==0:
            raise Exception('no matching first element')





    for ind, ExpTime in enumerate(ExpTimes):
        list_img= list(filter(lambda x: (x['TEXPMS']==ExpTime*1000 and x['ImgOrDark']==1),CCDitems_channel))
        list_dark=list(filter(lambda x: (x['TEXPMS']==ExpTime*1000 and x['ImgOrDark']==0),CCDitems_channel))
        
        
        #hej=CCDarray[np.where(filter(lambda x: (x['TEXPMS']==ExpTime), CCDarray))]
        #nparreven = nparr[np.where(nparr % 2 == 0)]
        
        #list_values = [ key for key,val in mydict.items() if val==value ]



        
        
        try: 
            CCDitem=list_img[0]
            CCDDitem1=list_dark[0]
            CCDDitem2=list_dark[1]
        except:
            raise Exception("list_img or list_dark is empty")
        

        image_no_dark=1.*CCDitem['IMAGE']-(0.5*CCDDitem2['IMAGE']+0.5*CCDDitem1['IMAGE'])
        

#        image_bias_sub = get_true_image(image_lsb, CCDitem.copy())
        image_desmeared = desmear_true_image(CCDitem.copy(),image_no_dark.copy())

#        image_dark_sub=subtract_dark(image_desmeared, CCDitem.copy())


        
        if ifile==0 : #No shutter
            ifileinfo='no shutter '
            CCDitem['image_sub']=image_no_dark/(ExpTime)
            CCDitem['image_processed']=image_desmeared/(ExpTime)
            

#            CCDplot(myfig,axel[ind],image_desmeared,'jrek')
  
            
        elif ifile==1: #Shutter
            ifileinfo='shutter '
            CCDitem['image_sub']=image_no_dark /(ExpTime-shuttertime)
            CCDitem['image_processed']=CCDitem['image_sub']#
        else:
            raise Exception('too many input files')
        list_scaled.append(CCDitem)
        
        #################################################
        # Plotting the results                          #
        #################################################
        
        
        sp=CCDplot(figs[0],axes[0][ind][ifile],CCDitem['image_sub'],title=ifileinfo+'Exptime= '+str(ExpTime))
        sp=CCDplot(figs[1],axes[1][ind][ifile],CCDitem['image_processed'],title=ifileinfo+'Exptime= '+str(ExpTime))

        if ifile==0:
            CCD_noshutter.append(CCDitem)
        elif ifile==1:
            CCD_shutter.append(CCDitem)
        else:
            raise Exception('too many input files')
    
     
        # if ifile==0: #plot dark cuirrent subtrr                    
        #     sp=CCDplot(figs[2],axes[2][ind][0],CCDitem['image_lsb'],title=ifileinfo+'Exp time = '+str(ExpTime))        
        #     sp=CCDplot(figs[2],axes[2][ind][1],image_dark_sub,title=ifileinfo+'Exp time = '+str(ExpTime))




for i in range(0,len(CCD_shutter)):
    ifile=2 #third row in plot
    
    diff_sub=(CCD_shutter[i]['image_sub']-CCD_noshutter[i]['image_sub'])      
    sp=CCDplot(figs[0],axes[0][i][ifile],diff_sub,title='(mid panel - midleft p.)')

    
    diff_lsb=(CCD_shutter[i]['image_processed']-CCD_noshutter[i]['image_processed'])       
    sp=CCDplot(figs[1],axes[1][i][ifile],diff_lsb,title='(mid panel - midleft p.)')
    

    ifile=3    #forth row in plot


    sp=plot_CCDimage_hmean(figs[0], axes[0][i][ifile],diff_sub[:,100:1800], title="horiz. mean midright p.")
    sp=plot_CCDimage_hmean(figs[1], axes[1][i][ifile],diff_lsb[:,100:1800], title="horiz. mean midright p.")
#    yax=range(0,diff_sub.shape[0])
#    sp=axes[0][i][ifile].plot(diff_sub[:,100:1800].mean(axis=1),yax)
#    axes[0][i][ifile].set_title("horiz. mean midright p.")
    
#    yax=range(0,diff_lsb.shape[0])
#    sp=axes[1][i][ifile].plot(diff_lsb[:,100:1800].mean(axis=1),yax)
#    axes[1][i][ifile].set_title("horiz. mean midright p.")


figs[0].savefig('test_desmearing'+str(SigM)+'_'+channel+'_noncompensated.jpg')
figs[1].savefig('test_desmearing'+str(SigM)+'_'+channel+'_compensated.jpg')
