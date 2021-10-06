#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 07:53:04 2020

@author: lindamegner
"""

import pandas as pd
import numpy as np


def readjorgfile(filename):

    df = pd.read_csv(filename, sep=" ", skiprows=(0), header=(0))
    
    df.columns=['ID','CCDout','ExpTime'] 
    
    df['ImgOrDark']=np.ones(df.shape[0])
    df['Measurement']=np.ones(df.shape[0])
    df.loc[::2, 'ImgOrDark']=0 # Set every second line to darkpic 
    
    #Remove data that is not measurements
    for index, row in df.iterrows():
        if len(row['ID'])<6 :
            df.iloc[[index],[4]]=0
            
          
    df_meas=df[df.Measurement==1]
    df_meas=df_meas.drop(columns=['Measurement', 'CCDout', 'ExpTime'])
    return df_meas