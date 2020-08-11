#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 08:57:37 2020

@author: lindamegner

Functions used to read in MATS images and data in different ways: From KTH, from Immage viewer and from rac files. 

"""

import numpy as np
import time
import matplotlib.pyplot as plt
import json
from PIL import Image

# some_file.py
import sys
# insert at 1, 0 is the script path (or '' in REPL)
#sys.path.insert(1, '/Users/lindamegner/MATS/retrieval/Level0/MATS-L0-processing-master')

import imagereader as imagereader
def readimageviewpic(dirname,picnr,rawflag):
# This function should not needed anymore.          
    if rawflag==1:
        imagefile= dirname +'rawoutput'+str(picnr) +'.pnm'
    else:
        imagefile= dirname +'output'+str(picnr) +'.pnm'
    txtfile=dirname +'output'+str(picnr) +'.txt'
    image_raw = np.float64(Image.open(imagefile))
    CCDitem=read_txtfile_create_CCDitem(txtfile)
    CCDitem['IMAGE']=image_raw    
    return CCDitem


def read_CCDitem_from_imgview(dirname,IDstring):
# This function used to be called readimageviewpic2    
# For many purposes it is probably better to use the fucntion read_pnm_image_and_txt      
    imagefile= dirname + '/'+ IDstring +'.pnm'
    txtfile= dirname +'/'+ IDstring +'_output.txt'
    try:
        image_raw = np.float64(Image.open(imagefile))
        CCDitem=read_txtfile_create_CCDitem(txtfile)
        CCDitem['IMAGE']=image_raw    
    except:
        print('There is something wrong with image file ',imagefile)
        CCDitem=-999
        raise Exception()

    return CCDitem
    
def read_pnm_image_and_txt(dirname,picid):
# Reads data from output image file and combines it with the txt file.
    imagefile= dirname +picid +'.pnm'
    txtfile=dirname +picid+'_output.txt'
    image_raw = np.float64(Image.open(imagefile))
    CCDitem=read_txtfile_create_CCDitem(txtfile)
    CCDitem['IMAGE']=image_raw    
    return CCDitem
  

def readimageviewpics(dirname,rawflag=0, filelist=[]):
    from os import listdir
    from os.path import isfile, join    
    #Reads  all images in directory
    all_files = [f for f in listdir(dirname) if isfile(join(dirname, f))]
    pnm_files = list(filter(lambda x: x[-4:] == '.pnm', all_files))
    CCDitems=[]
    for pnm_file in pnm_files:
        IDstring=pnm_file[:-4]
        CCDitem=read_CCDitem_from_imgview(dirname,IDstring,rawflag)
        if CCDitem != -999:
            CCDitems.append(CCDitem)
        
    return CCDitems


def readselectedimageviewpics(dirname,IDlist):
    CCDitems=[]
    for IDstring in IDlist:
        CCDitem=read_CCDitem_from_imgview(dirname,IDstring)
        if CCDitem != -999:
            CCDitems.append(CCDitem)
    return CCDitems






def read_MATS_image(filename, pathdir=''):
    # If you use an newer version of the rac extract reader then pthdir is not needed. Nwwer is approximately from the start of 2020. 
    json_file = open(filename,'r')
    CCD_image_data = json.load(json_file)
    json_file.close
    
    
    for item in CCD_image_data:
#        print(pathdir+str(CCD_image_data[i]['IMAGEFILE']) + '_data.npy')
        try:
            item['IMAGE'] = np.load(pathdir+str(item['imagefile']) + '_data.npy')
        except:    
            print('Warning, one image file seems corrupt and has been rmeoved')
            CCD_image_data.remove(item)
    
    return CCD_image_data


def read_CCDitem(rac_image_json_file,rac_sub_dir,itemnumber,rac_image_json_dir=''):
# reads data from one image (itemnumber) in the rac file

    CCD_image_data=read_MATS_image(rac_image_json_dir+rac_sub_dir+rac_image_json_file,rac_image_json_dir)

    CCDitem=CCD_image_data[itemnumber]
 

# The variables below are remain question marks
#    CCDitem['Noverflow'] = Noverflow # FBINOV?
 #   CCDitem['Ending'] = Ending       #not found in rac 
    


    
#    CCDitem['BC']=CCDitem['BC'][0] Bad culums ska vara en array
    CCDitem['CCDSEL']=CCDitem['CCDSEL'][0]
    CCDitem['EXPTS']=CCDitem['EXPTS'][0]
    CCDitem['EXPTSS']=CCDitem['EXPTSS'][0]
    CCDitem['FBINOV']=CCDitem['FBINOV'][0]
    CCDitem['FRAME']=CCDitem['FRAME'][0]
    CCDitem['JPEGQ']=CCDitem['JPEGQ'][0]
    CCDitem['LBLNK']=CCDitem['LBLNK'][0]
    CCDitem['NBC']=CCDitem['NBC'][0]
    CCDitem['NCOL']=CCDitem['NCOL'][0]
    CCDitem['NCSKIP']=CCDitem['NCSKIP'][0]
    CCDitem['NFLUSH']=CCDitem['NFLUSH'][0]
    CCDitem['NRBIN']=CCDitem['NRBIN'][0]        
    CCDitem['NROW']=CCDitem['NROW'][0]
    CCDitem['NRSKIP']=CCDitem['NRSKIP'][0]
    CCDitem['SID_mnemonic']=CCDitem['SID_mnemonic'][0]
    CCDitem['TBLNK']=CCDitem['TBLNK'][0]
    CCDitem['TEMP']=CCDitem['TEMP'][0]
    CCDitem['TEXPMS']=CCDitem['TEXPMS'][0]
    CCDitem['TIMING1']=CCDitem['TIMING1'][0] # Named Reserved1 in Georgis code /LM 20191115
    CCDitem['TIMING2']=CCDitem['TIMING2'][0] # Named Reserved2 in Georgis code /LM 20191115
    CCDitem['TIMING3']=CCDitem['TIMING3'][0] # Named VersionDate in Georgis code /LM 20191115
    CCDitem['VERSION']=CCDitem['VERSION'][0]
    CCDitem['WDWOV']=CCDitem['WDWOV'][0]
    CCDitem['ZERO']=CCDitem['ZERO'][0] 

    CCDitem['NColBinFPGA']=CCDitem['NColBinFPGA'][0]
    CCDitem['NColBinCCD']=CCDitem['NColBinCCD'][0]    
    CCDitem['DigGain']=CCDitem['DigGain'][0]
    CCDitem['TimingFlag']=CCDitem['TimingFlag'][0] 
    CCDitem['SigMode']=CCDitem['SigMode'][0]
    CCDitem['WinModeFlag']=CCDitem['WinModeFlag'][0]     
    CCDitem['WinMode']=CCDitem['WinMode'][0]  
    


    
    if int(CCDitem['CCDSEL'])==1:
        channel='IR1'
    elif int(CCDitem['CCDSEL'])==4:
        channel='IR2'
    elif int(CCDitem['CCDSEL'])==3:
        channel='IR3'
    elif int(CCDitem['CCDSEL'])==2:
        channel='IR4'
    elif int(CCDitem['CCDSEL'])==5:
        channel='UV1'
    elif int(CCDitem['CCDSEL'])==6:
        channel='UV2'
    elif int(CCDitem['CCDSEL'])==7:
        channel='NADIR'
    else:
        print('Error in CCDSEL, CCDSEL=',int(CCDitem['CCDSEL']))  
        
    
    CCDitem['channel']=channel
    


    return CCDitem


def read_CCDitems(rac_dir):
# reads data from one image (itemnumber) in the rac file

    CCD_image_data=read_MATS_image(rac_dir+'images.json')


    for CCDitem in CCD_image_data:
#        CCDitem=CCD_image_data[itemnumber]
     
    
    # The variables below are remain question marks
    #    CCDitem['Noverflow'] = Noverflow # FBINOV?
      #   CCDitem['Ending'] = Ending       #not found in rac 
        

    
        
    #    CCDitem['BC']=CCDitem['BC'][0] Bad culums ska vara en array
        CCDitem['CCDSEL']=CCDitem['CCDSEL'][0]
        CCDitem['EXPTS']=CCDitem['EXPTS'][0]
        CCDitem['EXPTSS']=CCDitem['EXPTSS'][0]
        CCDitem['FBINOV']=CCDitem['FBINOV'][0]
        CCDitem['FRAME']=CCDitem['FRAME'][0]
        CCDitem['JPEGQ']=CCDitem['JPEGQ'][0]
        CCDitem['LBLNK']=CCDitem['LBLNK'][0]
        CCDitem['NBC']=CCDitem['NBC'][0]
        CCDitem['NCOL']=CCDitem['NCOL'][0]
        CCDitem['NCSKIP']=CCDitem['NCSKIP'][0]
        CCDitem['NFLUSH']=CCDitem['NFLUSH'][0]
        CCDitem['NRBIN']=CCDitem['NRBIN'][0]        
        CCDitem['NROW']=CCDitem['NROW'][0]
        CCDitem['NRSKIP']=CCDitem['NRSKIP'][0]
        CCDitem['SID_mnemonic']=CCDitem['SID_mnemonic'][0]
        CCDitem['TBLNK']=CCDitem['TBLNK'][0]
        CCDitem['TEMP']=CCDitem['TEMP'][0]
        CCDitem['TEXPMS']=CCDitem['TEXPMS'][0]
        CCDitem['TIMING1']=CCDitem['TIMING1'][0] # Named Reserved1 in Georgis code /LM 20191115
        CCDitem['TIMING2']=CCDitem['TIMING2'][0] # Named Reserved2 in Georgis code /LM 20191115
        CCDitem['TIMING3']=CCDitem['TIMING3'][0] # Named VersionDate in Georgis code /LM 20191115
        CCDitem['VERSION']=CCDitem['VERSION'][0]
        CCDitem['WDWOV']=CCDitem['WDWOV'][0]
        CCDitem['ZERO']=CCDitem['ZERO'][0] 
    
        CCDitem['NColBinFPGA']=CCDitem['NColBinFPGA'][0]
        CCDitem['NColBinCCD']=CCDitem['NColBinCCD'][0]    
        CCDitem['DigGain']=CCDitem['DigGain'][0]
        CCDitem['TimingFlag']=CCDitem['TimingFlag'][0] 
        CCDitem['SigMode']=CCDitem['SigMode'][0]
        CCDitem['WinModeFlag']=CCDitem['WinModeFlag'][0]     
        CCDitem['WinMode']=CCDitem['WinMode'][0]  
        
    
    
        
        if int(CCDitem['CCDSEL'])==1:
            channel='IR1'
        elif int(CCDitem['CCDSEL'])==4:
            channel='IR2'
        elif int(CCDitem['CCDSEL'])==3:
            channel='IR3'
        elif int(CCDitem['CCDSEL'])==2:
            channel='IR4'
        elif int(CCDitem['CCDSEL'])==5:
            channel='UV1'
        elif int(CCDitem['CCDSEL'])==6:
            channel='UV2'
        elif int(CCDitem['CCDSEL'])==7:
            channel='NADIR'
        else:
            print('Error in CCDSEL, CCDSEL=',int(CCDitem['CCDSEL']))  
            
        
        CCDitem['channel']=channel
        
    

    return CCD_image_data

def ismember(a, b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, None) for itm in a]  # None can be replaced by any other "not in b" value



def readimg(filename):

#    filename='/Users/lindamegner/MATS/retrieval/Level1/data/June2018TVAC/IMAGES/1245852147_data.npy'
    data_arr=np.fromfile(filename, dtype='uint16')#check endianess, uint16 instead of ubit16 seems to work

 
    #convert header to binary
    #header_bin = [bin(np.double(data_arr[i])) for i in range(0,12)]#check if at least 16 bits, check if indexing correct
    header_bin = np.asarray([bin(data_arr[i]) for i in range(0,12)])#this is a string 
    #change format of header_bin elements to be formatted like matlab char array
    for i in range(0,len(header_bin)):
        header_bin[i]=header_bin[i][2:].zfill(16)
    #read header
    Frame_count = int(header_bin[0],2)
    NRow = int(header_bin[1],2)
    NRowBinCCD = int(header_bin[2][10:16],2)
    NRowSkip = int(header_bin[3][8:16],2)
    NCol = int(header_bin[4],2)
    NColBinFPGA = int(header_bin[5][2:8],2)
    NColBinCCD = int(header_bin[5][8:16],2)
    NColSkip = int(header_bin[6][5:16],2)
    N_flush = int(header_bin[7],2)
    Texposure_MSB = int(header_bin[8],2)
    Texposure_LSB = int(header_bin[9],2)
    Gain = int(header_bin[10],2)
    SignalMode = Gain >> 12 & 1   
    Temperature_read = int(header_bin[11],2)




    #read image
    if len(data_arr)< NRow*(NCol+1)/(2**(NColBinFPGA)):#check for differences in python 2 and 3
        img_flag = 0
        image = 0
        Noverflow = 0
        BlankLeadingValue = 0
        BlankTrailingValue = 0
        ZeroLevel = 0
        
        Reserved1 = 0
        Reserved2 = 0
        
        Version = 0
        VersionDate = 0
        NBadCol = 0
        BadCol = 0
        Ending = 'Wrong size'
    else:
        img_flag = 1
        image = np.reshape(np.double(data_arr[11+1:NRow*(NCol+1)+12]), (NCol+1,NRow))
        image = np.matrix(image).getH()
    
        #Trailer
        trailer_bin = np.asarray([bin(data_arr[i]) for i in range(NRow*(NCol+1)+12,len(data_arr))])
        Noverflow = int(trailer_bin[0],2)
        BlankLeadingValue = int(trailer_bin[1],2)
        BlankTrailingValue = int(trailer_bin[2],2)
        ZeroLevel = int(trailer_bin[3],2)
        
        Reserved1 = int(trailer_bin[4],2)
        Reserved2 = int(trailer_bin[5],2)
        
        Version = int(trailer_bin[6],2)
        VersionDate = int(trailer_bin[7],2)
        NBadCol = int(trailer_bin[8],2)
        BadCol = []
        Ending = int(trailer_bin[-1],2)
        
        if NBadCol > 0:
            BadCol = np.zeros(NBadCol)
            for k_bc in range(0,NBadCol):
                BadCol[k_bc] = int(trailer_bin[9+k_bc],2)#check if this actually works
    
    #original matlab code uses structured array, as of 20-03-2019 implementation as dictionary seems to be more useful choice
    #decision might depend on further (computational) use of data, which is so far unknown to me
    header = {}
    header['Size'] = len(data_arr)
    header['FRAME'] = Frame_count
    header['NROW'] = NRow
    header['NRBIN'] = NRowBinCCD
    header['NRSKIP'] = NRowSkip
    header['NCOL'] = NCol
    header['NColBinFPGA'] = NColBinFPGA
    header['NColBinCCD'] = NColBinCCD
    header['NCSKIP'] = NColSkip
    header['NFLUSH'] = N_flush
    header['TEXPMS'] = Texposure_LSB + Texposure_MSB*2**16
    header['DigGain'] = Gain & 0b1111
    header['SigMode'] = SignalMode
    header['TEMP'] = Temperature_read
    header['Noverflow'] = Noverflow # FBINOV?
    header['LBLNK'] = BlankLeadingValue
    header['TBLNK'] = BlankTrailingValue
    header['ZERO'] = ZeroLevel
    header['TIMING1'] = Reserved1
    header['TIMING2'] = Reserved2
    header['VERSION'] = Version
    header['TIMING3'] = VersionDate
    header['NC'] = NBadCol
    header['BC'] = BadCol
    header['Ending'] = Ending        
    

    
    header['CCDSEL'] = 1 #LM hardcoded - should not be used   
    header['channel'] = 'KTH test channel' 
    return image, header, img_flag


def readimage_create_CCDitem(path, file_number): #reads file from georigis stuff LM20191113

    

    
    filename = '%sF_0%02d/D_0%04d' % (path, np.floor(file_number/100),file_number)

    
#    filename='/Users/lindamegner/MATS/retrieval/Level1/data/June2018TVAC/IMAGES/1245852147_data.npy'
    data_arr=np.fromfile(filename, dtype='uint16')#check endianess, uint16 instead of ubit16 seems to work

 
    #convert header to binary
    #header_bin = [bin(np.double(data_arr[i])) for i in range(0,12)]#check if at least 16 bits, check if indexing correct
    header_bin = np.asarray([bin(data_arr[i]) for i in range(0,12)])#this is a string 
    #change format of header_bin elements to be formatted like matlab char array
    for i in range(0,len(header_bin)):
        header_bin[i]=header_bin[i][2:].zfill(16)
    #read header
    Frame_count = int(header_bin[0],2)
    NRow = int(header_bin[1],2)
    NRowBinCCD = int(header_bin[2][10:16],2)
    NRowSkip = int(header_bin[3][8:16],2)
    NCol = int(header_bin[4],2)
    NColBinFPGA = int(header_bin[5][2:8],2)
    NColBinCCD = int(header_bin[5][8:16],2)
    NColSkip = int(header_bin[6][5:16],2)
    N_flush = int(header_bin[7],2)
    Texposure_MSB = int(header_bin[8],2)
    Texposure_LSB = int(header_bin[9],2)
    Gain = int(header_bin[10],2)
    SignalMode = Gain >> 12 & 1      
    Temperature_read = int(header_bin[11],2)


    #read image
    if len(data_arr)< NRow*(NCol+1)/(2**(NColBinFPGA)):#check for differences in python 2 and 3
        img_flag = 0
        image = 0
        Noverflow = 0
        BlankLeadingValue = 0
        BlankTrailingValue = 0
        ZeroLevel = 0
        
        Reserved1 = 0
        Reserved2 = 0
        
        Version = 0
        VersionDate = 0
        NBadCol = 0
        BadCol = 0
        Ending = 'Wrong size'
    else:
        img_flag = 1
        image = np.reshape(np.double(data_arr[11+1:NRow*(NCol+1)+12]), (NCol+1,NRow))
        image = np.matrix(image).getH()
    
        #Trailer
        trailer_bin = np.asarray([bin(data_arr[i]) for i in range(NRow*(NCol+1)+12,len(data_arr))])
        Noverflow = int(trailer_bin[0],2)
        BlankLeadingValue = int(trailer_bin[1],2)
        BlankTrailingValue = int(trailer_bin[2],2)
        ZeroLevel = int(trailer_bin[3],2)
        
        Reserved1 = int(trailer_bin[4],2)
        Reserved2 = int(trailer_bin[5],2)
        
        Version = int(trailer_bin[6],2)
        VersionDate = int(trailer_bin[7],2)
        NBadCol = int(trailer_bin[8],2)
        BadCol = []
        Ending = int(trailer_bin[-1],2)
        
        if NBadCol > 0:
            BadCol = np.zeros(NBadCol)
            for k_bc in range(0,NBadCol):
                BadCol[k_bc] = int(trailer_bin[9+k_bc],2)#check if this actually works
    
    #original matlab code uses structured array, as of 20-03-2019 implementation as dictionary seems to be more useful choice
    #decision might depend on further (computational) use of data, which is so far unknown to me
    CCDitem = {}
#    CCDitem['Size'] = len(data_arr)
    CCDitem['FRAME'] = Frame_count
    CCDitem['NROW'] = NRow
    CCDitem['NRBIN'] = NRowBinCCD
    CCDitem['NRSKIP'] = NRowSkip
    CCDitem['NCOL'] = NCol
    CCDitem['NColBinFPGA'] = NColBinFPGA 
    CCDitem['NColBinCCD'] = NColBinCCD
    CCDitem['NCSKIP'] = NColSkip
    CCDitem['NFLUSH'] = N_flush
    CCDitem['TEXPMS'] = Texposure_LSB + Texposure_MSB*2**16
    CCDitem['DigGain'] = Gain & 0b1111
    CCDitem['SigMode'] = SignalMode 
    CCDitem['TEMP'] = Temperature_read
    CCDitem['Noverflow'] = Noverflow # FBINOV?
    CCDitem['LBLNK'] = BlankLeadingValue  
    CCDitem['TBLNK'] = BlankTrailingValue
    CCDitem['ZERO'] = ZeroLevel
    CCDitem['TIMING1'] = Reserved1
    CCDitem['TIMING2'] = Reserved2 
    CCDitem['VERSION'] = Version
    CCDitem['TIMING3'] = VersionDate 
    CCDitem['NBC'] = NBadCol
    CCDitem['BC'] = BadCol
    CCDitem['Ending'] = Ending       #not found in rac 
    
    CCDitem['IMAGE']=image
    
    CCDitem['CCDSEL'] = 1 # Note that this is incorrect but hte KTH test CCD is unknown    
    CCDitem['channel'] = 'KTH test channel' 
    

    
  


    
    return CCDitem, img_flag



# Why has image values of 10 000 when raw_image has values of 400?
def readracimg(filename):
#   Linda Megner; function to read in from rac file but yield similar result 
#    as when read in by readimg . Note that this header has more info.
    
    
    
    image, metadata = imagereader.read_MATS_image(filename) 
    image=np.float64(image)
#   LM Do we want image to be float or integer? When do we want to convert it?     

    header = {}
#    header['Size'] = image.size
#   LM: Note this is NOT equivalent to what is in readim    
#    header['Size'] = len(data_arr). This variable nopt needed

  
    header['Frame_count'] = metadata['FRAME'][0]
    header['NRow'] = metadata['NROW'][0]
    header['NRowBinCCD'] = metadata['NRBIN'][0]
    header['NRowSkip'] = metadata['NRSKIP'][0]
    header['NCol'] = metadata['NCOL'][0]
    header['NColBinFPGA'] = metadata['NCBIN'][0] >> 8 & 0b1111
    header['NColBinCCD'] = metadata['NCBIN'][0] & 0b11111111
    header['NColSkip'] = metadata['NCSKIP'][0]
    header['N_flush'] = metadata['NFLUSH'][0]
    header['Texposure'] = metadata['TEXPMS'][0]
    header['DigGain'] = metadata['GAIN'] & 0b1111
    header['SignalMode'] = metadata['GAIN'] >> 12 & 1          
     
    
    




    
    header['Temperature_read'] = metadata['TEMP'][0] #LM what is this the temperature of and how do I convert it? 
#  Temperature is measured in the CRBD – readout ADC, so closer to OBC.
#It’s a 15 bit value from 2.048V (+/-0.5%) reference.
#So 4307/32768*2048 = 269 mV
#According to datasheet ADC outputs 273mV @ 25°C with 0.85 mV/°C sensitivity. So it’s approximately 20-21°C when the image was taken.
    header['Noverflow'] =  metadata['FBINOV'][0]  
    header['BlankLeadingValue'] =  metadata['LBLNK'][0]  
    header['BlankTrailingValue'] =  metadata['TBLNK'][0]  
    header['ZeroLevel'] =  metadata['ZERO'][0] 

    
    header['Reserved1'] = metadata['TIMING1'][0] #Check that this is correct by comparing old and new software ICD
    header['Reserved2'] = metadata['TIMING2'][0] #Check that this is correct by comparing old and new software ICD
    header['Version'] = metadata['VERSION'][0]
    header['VersionDate'] = metadata['TIMING3'][0] #Check that this is correct by comparing old and new software ICD
    header['NBadCol'] = metadata['NBC'][0]
    header['BadCol'] = metadata['BC']
#    header['Ending'] = 64175    NO INFO IN RAC FILE metadata 
    
#   More information not used in readimg
    header['CCDSEL'] = metadata['CCDSEL'][0]  
    header['EXPTS'] = metadata['EXPTS'][0]   
    header['EXPTSS'] = metadata['EXPTSS'][0]   
    header['JPEGQ'] = metadata['JPEGQ'][0]  
    header['SID_mnemonic'] = metadata['SID_mnemonic'][0] 
    header['WDW'] = metadata['WDW'][0]   
    header['WDWOV'] = metadata['WDWOV'][0]
    
    header['TimingFlag'] =  metadata['GAIN'] >> 8 & 1
  
    
#    img_flag=1 #LM is this needed? Ask Georgi
    return image, header





def readimgpath(path, file_number, plot):
    #stuff happens

    filename = '%sF_0%02d/D_0%04d' % (path, np.floor(file_number/100),file_number)

    image, header, img_flag = readimg(filename)



    if plot>0:
        #do the plotting
        mean_img=np.mean(image)
    
        plt.imshow(image, cmap='viridis', vmin=mean_img-100, vmax=mean_img+100)
        plt.title('CCD image')
        plt.xlabel('Pixels')
        plt.ylabel('Pixels')
        plt.show()#equivalent to hold off
        
        print(header)
        time.sleep(0.1)
    
    return image, header


def read_txtfile(filepath):    
    dict = {}
    with open(filepath, "r") as f:
        for line in f:
            (key, val) = line.split("=")
            if key=='id':
                dict['ID'] = val[1:-2] #For some reason a space has been added before and an ewnd og line character in the end - remove these
            else:
                dict[key] = int(val)
    return dict
    
def read_txtfile_create_CCDitem(filepath):
    # function used to read txt output from image viewer
    CCDitem=read_txtfile(filepath)


# Extract variables from certain bits within the same element, see 6.4.1 Software ICD /LM 20191115               

    CCDitem['NColBinFPGA'] =  CCDitem['NCBIN'] >> 8 & 0b1111
    CCDitem['NColBinCCD'] = CCDitem['NCBIN'] & 0b11111111
    del CCDitem['NCBIN']
    CCDitem['DigGain'] = CCDitem['GAIN'] & 0b1111
    CCDitem['TimingFlag'] = CCDitem['GAIN'] >> 8 & 1
    CCDitem['SigMode'] =CCDitem['GAIN'] >> 12 & 1          
    del CCDitem['GAIN']
    CCDitem['WinModeFlag']=CCDitem['WDW'] >> 7 & 1
    CCDitem['WinMode']=CCDitem['WDW'] & 0b111
    del CCDitem['WDW']
    # Hack to make image viewing output the same as the rac file outputs
    CCDitem['id']=CCDitem['ID']
    del CCDitem['ID']
    
    
    print('The following values may be incorrect:  NColBinFPGA, NColBinCCD DigGain TimingFlag SigMode inModeFlag WinMode')
    if int(CCDitem['CCDSEL'])==1:
        channel='IR1'
    elif int(CCDitem['CCDSEL'])==4:
        channel='IR2'
    elif int(CCDitem['CCDSEL'])==3:
        channel='IR3'
    elif int(CCDitem['CCDSEL'])==2:
        channel='IR4'
    elif int(CCDitem['CCDSEL'])==5:
        channel='UV1'
    elif int(CCDitem['CCDSEL'])==6:
        channel='UV2'
    elif int(CCDitem['CCDSEL'])==7:
        channel='NADIR'
    else:
        print('Error in CCDSEL, CCDSEL=',int(CCDitem['CCDSEL']))  
    CCDitem['channel']=channel
    
    return CCDitem