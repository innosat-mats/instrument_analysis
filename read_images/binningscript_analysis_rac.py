# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:23:19 2020

@author: bjorn
"""


# fmt: off
import sys
sys.path.insert(0, "/home/olemar/Projects/MATS/MATS-L0-processsing")
sys.path.insert(1, "/home/olemar/Projects/MATS/MATS-L1-processsing")


import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from LindasCalibrationFunctions import plotCCDitem
from LindasCalibrationFunctions import plot_simple
from L1_calibration_functions import desmear_true_image
import copy
import read_in_functions
from read_in_functions import read_CCDitem_from_imgview, readimageviewpics, read_MATS_image

# fmt: on


# TO DO: COMBINE THE FOLLOWING TWO FUNCITONS INTO ONE THAT BINS ACCORDING TO
# BOTH FPGA AND ON-CHIP & ROW SETTINGS;

def bin_ref(ref, ccd):

    # simple code for binning

    nrow, ncol, nrskip, ncskip, nrbin, ncbin, exptime = (ccd['NROW'], ccd['NCOL']+1,
                                                         ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NColBinCCD'], ccd['TEXPMS'])

    nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr, exptimer = (ref['NROW'], ref['NCOL']+1,
                                                                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NColBinCCD'], ref['TEXPMS'])

    exptimefactor = int((exptime-2000)/(exptimer-2000))
    # reference image that will be binned according to 'ccd' settings
    imgref = ref['IMAGE']

    # in case reference image is already a binned image
    ncbin, nrbin = int(ncbin/ncbinr), int(nrbin/nrbinr)

    # images must cover the same ccd section
    if ncskip == ncskipr and nrskip == nrskipr:

        colbin = np.zeros([nrowr, ncol])

        for j in range(0, ncol):
            colbin[:, j] = imgref[:, j*ncbin:(j+1)*ncbin].sum(axis=1)

        # declare zero array for row binning
        binned = np.zeros([nrow, ncol])

        for j in range(0, nrow):
            binned[j, :] = colbin[j*nrbin:(j+1)*nrbin, :].sum(axis=0)

        binned = binned*exptimefactor
        return binned

    else:

        sys.exit('Error: images not from the same CCD region.')


def bin_ref_FPGA(ref, ccd):

    # simple code for binning
    nrow, ncol, nrskip, ncskip, nrbin, ncbin, exptime = (ccd['NROW'], ccd['NCOL']+1,
                                                         ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NColBinCCD'], ccd['TEXPMS'])

    nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr, exptimer = (ref['NROW'], ref['NCOL']+1,
                                                                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NColBinCCD'], ref['TEXPMS'])

    exptimefactor = int((exptime-2000)/(exptimer-2000))
    # reference image that will be binned according to 'ccd' settings
    imgref = ref['IMAGE']

    # in case reference image is already a binned image
    ncbin, nrbin = int(ncbin/ncbinr), int(nrbin/nrbinr)

    # images must cover the same ccd section
    if ncskip == ncskipr and nrskip == nrskipr:

        colbin = np.zeros([nrowr, ncol])

        for j in range(0, ncol):
            colbin[:, j] = imgref[:, j*ncbin:(j+1)*ncbin].sum(axis=1)

        # declare zero array for row binning
        binned = np.zeros([nrow, ncol])

        for j in range(0, nrow):
            binned[j, :] = colbin[j*nrbin:(j+1)*nrbin, :].sum(axis=0)

        binned = binned*exptimefactor
        return binned

    else:

        sys.exit('Error: images not from the same CCD region.')


def img_diff(image1, image2):

    return image1-image2


####################################################
#############        LOAD DATA      ################
####################################################

def main():  # type of binning

    for channel in range(1, 8):

        dirname = (
            '/home/olemar/Projects/MATS/MATS-data/binning_test_20200812_racfiles/binning/')
        CCDitems = []
        IDstrings = []
        binned = []

        os.chdir(dirname)

        CCDitems = read_in_functions.read_CCDitems(
            dirname
        )

        CCDitems_use = []

        for i in range(len(CCDitems)):
            if CCDitems[i]['CCDSEL'] == channel:
                CCDitems_use.append(CCDitems[i])

        CCDitems = CCDitems_use
        # long exposure, short exposure, new reference images
        CCDl_list = np.copy(CCDitems[0::4])
        CCDs_list = np.copy(CCDitems[1::4])
        CCDr_list = np.copy(CCDitems[2::4])
        CCDrs_list = np.copy(CCDitems[3::4])  # reference short (never binned)

    # SUBTRACT DARK IMAGES

        CCDl_sub_img, CCDr_sub_img = [], []
        test_type = np.array([])
        for i in range(0, len(CCDs_list)):

            print(i)
            if CCDl_list[i]['TEXPMS'] != CCDr_list[0]['TEXPMS']:
                test_type = np.append(test_type, 'exp')
            elif CCDl_list[i]['NCBIN CCDColumns'] != CCDr_list[0]['NCBIN CCDColumns']:
                test_type = np.append(test_type, 'col')
            elif CCDl_list[i]['NRBIN'] != CCDr_list[0]['NRBIN']:
                test_type = np.append(test_type, 'row')
            else:
                test_type = np.append(test_type, 'ref')
        #    subtract dark current from both long and references
            CCDl_sub_img.append(
                img_diff(CCDl_list[i]['IMAGE'].copy(), CCDs_list[i]['IMAGE'].copy()))
            CCDr_sub_img.append(img_diff(
                CCDr_list[i]['IMAGE'].copy(), CCDrs_list[i]['IMAGE'].copy()))  # update

        # # plot the images with removed dark current
        # fig1, axs1 = plt.subplots(1, len(CCDs_list), figsize=(
        #     15, 6), facecolor='w', edgecolor='k')
        # fig1.subplots_adjust(hspace=1, wspace=.001)
        # axs1 = axs1.ravel()
        # fig1.suptitle('images with subtracted dark images')

        # for i in range(0, len(CCDs_list)):
        #     im = axs1[i].imshow(CCDl_sub_img[i], cmap='jet')
        #     mean = CCDl_sub_img[i].mean()
        #     std = CCDl_sub_img[i].std()
        #     fig1.colorbar(im, ax=axs1[i])
        #     im.set_clim(mean-2*std, mean+2*std)
        #     axs1[i].set_title('long - short: ' + str(CCDs_list[i]['NRBIN']) + 'x'
        #                       + str(CCDs_list[i]['NColBinCCD']*2**CCDs_list[i]['NColBinFPGA']))

        ####################################################
        #############        BINNING       #################
        ####################################################

        # copy settings from long image (for binning settings)
        bin_input = copy.deepcopy(CCDl_list)

        # replace images with the images with subtrated dark (I think this is old and not used)
        for i in range(0, len(CCDs_list)):
            bin_input[i]['IMAGE'] = CCDl_sub_img[i].copy()

        # create manually binned images
        for i in range(0, len(CCDs_list)):

            # update reference image that should be binned manually
            ref = copy.deepcopy(CCDr_list[i])
            ref['IMAGE'] = CCDr_sub_img[i].copy()

            # bin reference image according to bin_input settings
            binned.append(bin_ref(copy.deepcopy(ref), bin_input[i].copy()))

        # # plot histograms of manual and instrument bins
        # fig3, axs3 = plt.subplots(3, len(CCDs_list), figsize=(
        #     15, 6), facecolor='w', edgecolor='k')
        # fig3.subplots_adjust(hspace=1, wspace=.001)
        # axs3 = axs3.ravel()

        # # start and stop for plotting
        # rstart, rstop = 1, -1
        # cstart, cstop = 1, -1

        # # bin (histogram) size
        # binn = 80

        # for i in range(0, len(CCDs_list)):

        #     print(i)

        #     inst_bin = CCDl_sub_img[i].copy()

        #     mean_man = binned[i].mean()
        #     std_man = binned[i].std()

        #     mean_inst = inst_bin.mean()
        #     std_inst = inst_bin.std()

        #     axs3[i].hist(binned[i][rstart:rstop, cstart:cstop].ravel(), range=(mean_man-3*std_man, mean_man+3*std_man), bins=binn, alpha=0.6,
        #                  color="skyblue", label='manual', density=True)
        #     axs3[i+len(CCDs_list)].hist(inst_bin[rstart:rstop, cstart:cstop].ravel(), range=(mean_inst-3*std_man, mean_inst+3*std_inst), bins=binn, alpha=0.4,
        #                                 color='red', label='instrument', density=True)
        #     axs3[i+2*len(CCDs_list)].hist(inst_bin[rstart:rstop, cstart:cstop].ravel(), bins=binn, range=(mean_man-3*std_man, mean_man+3*std_man), alpha=0.4,
        #                                   color='red', label='instrument', density=True)
        #     axs3[i+2*len(CCDs_list)].hist(binned[i][rstart:rstop, cstart:cstop].ravel(), bins=binn, range=(mean_man-3*std_man, mean_man+3*std_man), alpha=0.6,
        #                                   color='skyblue', label='manual', density=True)
        #     axs3[i].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
        #                       + str(CCDs_list[i]['NColBinCCD']*2**CCDs_list[i]['NColBinFPGA'])+' man')
        #     axs3[i+len(CCDs_list)].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
        #                                      + str(CCDs_list[i]['NColBinCCD']*2**CCDs_list[i]['NColBinFPGA']) + ' inst')
        #     axs3[i+2*len(CCDs_list)].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
        #                                        + str(CCDs_list[i]['NColBinCCD']*2**CCDs_list[i]['NColBinFPGA']) + ' comp')

        mean_man_tot = np.array([])
        mean_inst_tot = np.array([])

        for i in range(0, len(CCDs_list)):
            inst_bin = CCDl_sub_img[i].copy()
            mean_man_tot = np.append(mean_man_tot, binned[i].mean())
            mean_inst_tot = np.append(mean_inst_tot, inst_bin.mean())

        plt.plot(mean_man_tot[test_type == 'exp'],
                 mean_inst_tot[test_type ==
                               'exp'], '.', mean_man_tot[test_type == 'row'],
                 mean_inst_tot[test_type ==
                               'row'], '+', mean_man_tot[test_type == 'col'],
                 mean_inst_tot[test_type == 'col'], 'x')
    plt.ylabel('measured values')
    plt.xlabel('simulated values')
    plt.show()


if __name__ == "__main__":
    # execute only if run as a script
    main()
