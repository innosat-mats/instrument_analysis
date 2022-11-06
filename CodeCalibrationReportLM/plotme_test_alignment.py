#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 22:52:21 2022

@author: lindamegner

"""
laser1_start=54
laser1_stop=59
laser2_start=41
laser2_stop=49
deuterium_LFOV_start=40
deuterium_LFOV_stop=46
deuterium_RFOV_start=47
deuterium_RFOV_stop=53


for ind, CCDitem in enumerate(CCDitems[laser2_start:laser2_stop]):
    (
        image_lsb,
        image_bias_sub,
        image_desmeared,
        image_dark_sub,
        image_calib_nonflipped,
        image_calibrated,
        image_common_fov,
        errors
    ) = L1_calibrate(CCDitem, instrument)

    if plot:
        fig, ax = plt.subplots(6, 1)
        plot_CCDimage(image_lsb, fig, ax[0], "Original LSB")
        plot_CCDimage(image_bias_sub, fig, ax[1], "Bias subtracted")
        plot_CCDimage(image_desmeared, fig, ax[2], " Desmeared LSB")
        sp=plot_CCDimage(image_dark_sub, fig, ax[3], " Dark current subtracted LSB")
        fix_clim=sp.get_clim()
        plot_CCDimage(image_calib_nonflipped, fig, ax[4], " Flat field compensated LSB", clim=fix_clim)
        plot_CCDimage(image_common_fov, fig, ax[5], " Image common field of view", clim=fix_clim)
        fig.suptitle(CCDitem["channel"])
            
    plot_CCDimage(image_lsb, fig1, ax1[ind], CCDitem['channel'])    
    fig.savefig('images/Alignment_'+CCDitem['channel']+'.jpg')             
    sp=plot_CCDimage_transp(image_common_fov, fig2, ax2, 'All channels', clim=fix_clim, alpha=0.3)    
    
    