#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:53:36 2022

@author: lindamegner


This is more or less a duplicate of read_and_calibrate_all_files_in_directory but it is meant to be used as a script.
"""


import matplotlib.pyplot as plt
from mats_l1_processing.experimental_utils import plot_CCDimage



    
def select_CCDitems(CCDitems, key, value): 
    """
    Function to seache through image items to find specific settings

    Parameters
    ----------
    CCDitems : list
        DESCRIPTION.
    key : string: name og CCDsetting or quanitity, i.e. CCDitem dictonary intex, e.g. 'channel'
        DESCRIPTION.
    value : Value of CCDsetting or quanitity, e.g. 'IR4'
        DESCRIPTION.

    Returns
    -------
    CCDitems_select: List of CCDitems

    """

    CCDitems_select= list(filter(lambda x: ( x[key]==value),CCDitems))
    return CCDitems_select

def select_CCDitems_using_list(CCDitems, key, valuelist):

    
    for value in valuelist:
        CCDitems=select_CCDitems(CCDitems, key, value)
    
    return CCDitems



def select_CCDitems_using_keyvaluedict(CCDitems, key_value_dict):
    """
    
    Parameters
    ----------
    CCDitems : LIST of CCDitems
    key_value_dict : DICTIONARY specifying the keys to sort on and their values

    Returns
    -------
    LIST if CCDitems that pass all selection criteria in key_value_list

    """
    CCDitems_select=CCDitems
    for key, value in key_value_dict.items():
        CCDitems_select=select_CCDitems(CCDitems_select,key,value)
    
    return CCDitems_select



def sort_images_in_dirs(CCDitems, key_value_dict, clim=999, path='.', whattoplot="IMAGE"):
    """
    Takes a CCDitems list, plots the images and sorts the plots in directories according to the key_value_dict
    

    Parameters
    ----------
    CCDitems : LIST of CCDitems
    key_value_dict : DICTIONARY specifying the keys to sort on and their values
    clim : optional
        DESCRIPTION. The default is "999".
    path : STRING, optional
        DESCRIPTION. The default is ".".
    whattoplot : STRING, optional
        DESCRIPTION. The default is "IMAGE".

    Returns
    -------
    None.

    """
    import os
    
    for key, value in key_value_dict.items():
        keyval=str(key)+'_'+str(value)
        dirname=path+'/'+keyval
        CCDitems_select=select_CCDitems(CCDitems,key,value)
        os.mkdir(dirname)
        plot_CCDitems(CCDitems_select, path=dirname, title=keyval, whattoplot=whattoplot)
        
 
        
 
def sort_images_plot(CCDitems, key_value_dict, clim=999, path='.', whattoplot="IMAGE"):
    """
    Takes a CCDitems list, and plots selected images according to the key_value_dict
    

    Parameters
    ----------
    CCDitems : LIST of CCDitems
    key_value_dict : DICTIONARY specifying the keys to sort on and their values
    clim : optional
        DESCRIPTION. The default is "999".
    path : STRING, optional
        DESCRIPTION. The default is ".".
    whattoplot : STRING, optional
        DESCRIPTION. The default is "IMAGE".    


    Returns
    -------
    None.

    """
    
    keyval=''
    for key, value in key_value_dict.items():

        keyval=keyval+'_'+str(key)+'_'+str(value)

        #dirname=path+'/'+keyval
        CCDitems=select_CCDitems(CCDitems,key,value)
        #os.mkdir(dirname)
        #plot_CCDitems(CCDitems_select, path=dirname, title=keyval)
    if len(CCDitems)>10:
        raise Exception('Too many imgaes to plot in one plot - make a more narrpw selection')
    nr_of_plots=len(CCDitems)
    if len(CCDitems)==1: nr_of_plots=1
    fig, ax= plt.subplots(nr_of_plots,1)
    for ind, CCDitem in enumerate(CCDitems):

        sp=plot_CCDimage(CCDitem[whattoplot], fig, ax[ind], CCDitem['id']+'_'+keyval,clim)
            #fig.suptitle(channel)
            #sp.set(visible='False')        
    fig.savefig(path+'/image_'+CCDitem['id']+'.jpg')
            #plt.close(fig)
        
        
def create_plot_directory_tree(CCDitems, key_value_dict, clim=999, path='.'):
    """
    NOT FUNCTIONING YET - this function may not be needed
    
    Takes a CCDitems list, plots the images and sorts the plots in directories according to the key_value_dict
    
    

    Parameters
    ----------
    CCDitems : LIST of CCDitems
    key_value_dict : DICTIONARY specifying the keys to sort on and their values
    path : STRING, optional
        DESCRIPTION. The default is ".".


    Returns
    -------
    None.

    """
    import os
    
    dirname=path
    for key, value in key_value_dict.items():

        keyval=str(key)+'_'+str(value)
    
        dirname=dirname+'/'+keyval
        os.mkdir(dirname)
        CCDitems=select_CCDitems(CCDitems,key,value)
        
        
        plot_CCDitems(CCDitems, path=dirname, title=keyval) 
        
        
        

def plot_CCDitems(CCDitems, title="", clim=999, aspect="auto", path=".", whattoplot="IMAGE"):
    """
    Plots all CCDitems in the directory given as input

    Parameters
    ----------
    CCDitems : LIST OF DICTS
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".
    clim : TYPE, optional
        DESCRIPTION. The default is 999.
    aspect : TYPE, optional
        DESCRIPTION. The default is "auto".
    path : STRING, optional
        DESCRIPTION. The default is ".".
    whattoplot : STRING, optional
        DESCRIPTION. The default is "IMAGE".

    Returns
    -------
    None.

    """

    for CCDitem in CCDitems:
        fig = plt.figure()
        ax = fig.gca()
        sp=plot_CCDimage(CCDitem[whattoplot], fig, ax, CCDitem['id']+'_'+title,clim,aspect)
        #fig.suptitle(channel)
        #sp.set(visible='False')        
        fig.savefig(path+'/image_'+CCDitem['id']+'.jpg')
        #plt.close(fig)

        
