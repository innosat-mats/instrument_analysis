#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""



from mats_l1_processing.experimental_utils import  plotCCDitem, read_all_files_in_protocol
from mats_l1_processing.experimental_utils import readprotocol


from mats_l1_processing.read_in_functions import read_all_files_in_root_directory
import matplotlib.pyplot as plt



directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/LaserAlignmentTest_RacFiles0909to0910Duplicated/'



protocol1='protocol_lasertest_morning20210909.txt'

protocol2='protocol_alignment_w_IR3_morning2021010.txt'

#directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/BinningFlatfieldsIR3_210910/'
#protocol='readprotocol_IR3.txt'

protocol=protocol1

read_from="rac" 

#fig3, axes3 = plt.subplots(2,1,figsize=(10,10))

readfromprot=True
if readfromprot:
    fig, ax=plt.subplots(5,1)
    df_protocol=readprotocol(directory+protocol)
    CCDitems=read_all_files_in_protocol(df_protocol, read_from,directory)
    if protocol==protocol2:
        CCDitems_run1=CCDitems[0:4]
        CCDitems_run2=CCDitems[4:8]
        CCDitems_run3=CCDitems[8:12]
        CCDitems_run4=CCDitems[12:16]
        CCDitems=CCDitems_run2

        

    for i, CCDitem in enumerate(CCDitems):

         plotCCDitem(CCDitem,fig, ax[i], title=CCDitem['channel']+ ' '+CCDitem['id'])

   # fig.savefig('images/lamp_alignment_IRchannels_210910_run2.jpg') 
    
else:

    CCDitems=read_all_files_in_root_directory(read_from,directory)
#calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

    CCDitems_short = []
    for CCDitem in CCDitems[250:270]: #250-270 appears to be some kind of laser points 
        CCDitems_short.append(CCDitem)
             

    for CCDitem in CCDitems_short:
         fig=plt.figure()
         ax=fig.gca()
         plotCCDitem(CCDitem,fig, ax, title=CCDitem['channel']+ ' '+CCDitem['id'])
         print(CCDitem['channel']+' '+CCDitem['id'])
    
        
    
    
