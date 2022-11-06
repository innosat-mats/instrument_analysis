#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 09:38:11 2022

@author: lindamegner
"""



def meanbin_image_with_BC(CCDitem, image_nonbinned=None):
    """
    This is a function to mean-bin an image (taking nskip into account) without any offset or blanks. Bad columns are skipped.
    This code is used for an image which has already been treated with regards to offset and BC(get_true_image)
    For instance when binning the flatfield image. Code is modified from Georgis predict_image 


    Args:
        CCDitem:  dictonary containing CCD image and information
        image_nonbinned (optional): numpy array image

    Returns: 
        meanbinned_image: binned image (by taking the average) according to the info in CCDitem 

    """

            
    if image_nonbinned is None:
        image_nonbinned = CCDitem["IMAGE"]
        
    #Check if image needs to be binned or shifted    
    totbin=int(CCDitem["NRBIN"])*int(CCDitem["NCBIN CCDColumns"])*int(CCDitem["NCBIN FPGAColumns"])  
    if (totbin>1 or CCDitem["NCSKIP"] >0 or CCDitem["NRSKIP"]>0):


        ncol = int(CCDitem["NCOL"]) + 1
        nrow = int(CCDitem["NROW"])
    
        nrowskip = int(CCDitem["NRSKIP"])
        ncolskip = int(CCDitem["NCSKIP"])
    
        nrowbin = int(CCDitem["NRBIN"])
        ncolbinC = int(CCDitem["NCBIN CCDColumns"])
        ncolbinF = int(CCDitem["NCBIN FPGAColumns"])
    
        bad_columns = CCDitem["BC"]
    
        if nrowbin == 0:  # no binning means beaning of one
            nrowbin = 1
    
        if ncolbinC == 0:  # no binning means beaning of one
            ncolbinC = 1
    
        if ncolbinF == 0:  # no binning means beaning of one
            ncolbinF = 1
    
        ncolbintotal = ncolbinC * ncolbinF
    
    
        image = np.zeros((nrow, ncol))  # no offset
        nr_of_entries = np.zeros((nrow, ncol))
    
        for j_r in range(0, nrow):  # check indexing again
            for j_c in range(0, ncol):
                for j_br in range(0, nrowbin):  # account for row binning on CCD
                    for j_bc in range(0, ncolbintotal):  # account for column binning
                        # LM201030: Go through all unbinned columns(both from FPGA and onchip) that belongs to one superpixel(j_r,j_c) and if the column is not Bad, add the signal of that unbinned pixel to the superpixel (j_r,j_c)
                        # out of reference image range
                        if (j_r) * nrowbin + j_br + nrowskip > 511:
                            break
                        elif (j_c) * ncolbinC * ncolbinF + j_bc + ncolskip > 2048:
                            break
    
                        if (
                            ncolbinC > 1
                            and (j_c) * ncolbinC * ncolbinF + j_bc + ncolskip in bad_columns
                        ):  # +1 becuase Ncol is +1
                            continue
                        else:
    
                            # add only the actual signal from every pixel 
                            image[j_r, j_c] = (
                                image[j_r, j_c]  
                                + image_nonbinned[
                                    (j_r) * nrowbin + j_br + nrowskip,
                                    (j_c) * ncolbinC * ncolbinF + j_bc + ncolskip,
                                ]  
                            )
    
                            nr_of_entries[j_r, j_c] = nr_of_entries[j_r, j_c] + 1
    
        meanbinned_image = image / nr_of_entries
    
    
    else:
        meanbinned_image=image_nonbinned

    
    return meanbinned_image

    

def bin_image_using_predict(header, reference_image="999"):
    """
    this is a function to predict an image read out from the CCD with a given set
    of parameters, based on a reference image (of size 511x2048)
    """

    ncol = int(header["NCOL"]) + 1
    nrow = int(header["NROW"])

    nrowskip = int(header["NRSKIP"])
    ncolskip = int(header["NCSKIP"])

    nrowbin = int(header["NRBIN"])
    ncolbinC = int(header["NCBIN CCDColumns"])
    ncolbinF = int(header["NCBIN FPGAColumns"])

    blank = int(header["TBLNK"])

    blank_off = blank - 128

    # gain=2**(int(header['Gain']) & 255) #use for old data format
    gain = 2.0 ** header["GAIN Truncation"]
    bad_columns = header["BC"]

    if nrowbin == 0:  # no binning means beaning of one
        nrowbin = 1

    if ncolbinC == 0:  # no binning means beaning of one
        ncolbinC = 1

    if ncolbinF == 0:  # no binning means beaning of one
        ncolbinF = 1

    ncolbintotal = ncolbinC * ncolbinF

    if type(reference_image) == str:
        if reference_image == "999":
            reference_image = header["IMAGE"]
        else:
            raise ValueError('Invalid reference image')

    # bad column analysis
    n_read, n_coadd = binning_bc(ncol, ncolskip, ncolbinF, ncolbinC, header["BC"])

    image = np.zeros((nrow, ncol))
    image[:, :] = 128  # offset

    finished_row = 0
    finished_col = 0
    for j_r in range(0, nrow):  # check indexing again
        for j_c in range(0, ncol):
            for j_br in range(0, nrowbin):  # account for row binning on CCD
                if j_br == 0:
                    image[j_r, j_c] = (
                        image[j_r, j_c] + n_read[j_c] * blank_off
                    )  # here we add the blank value, only once per binned row
                    # (LM201025 n_read is the number a superbin has been devided into to be read. So if no badcolums or fpga binning then n_read=1.
                for j_bc in range(0, ncolbintotal):  # account for column binning
                    # LM201030: Go through all unbinned columns(both from FPGA and onchip) that belongs to one superpixel(j_r,j_c) and if the column is not Bad, add the signal of that unbinned pixel to the superpixel (j_r,j_c)
                    # out of reference image range
                    if (j_r) * nrowbin + j_br + nrowskip > 511:
                        break
                    elif (j_c) * ncolbinC * ncolbinF + j_bc + ncolskip > 2048:
                        break

                    # removed +1 after bad_columns, unclear why it was added
                    # TODO
                    if (
                        ncolbinC > 1
                        and (j_c) * ncolbinC * ncolbinF + j_bc + ncolskip in bad_columns
                    ):  # +1 becuase Ncol is +1
                        continue
                    else:

                        # add only the actual signal from every pixel (minus blank)
                        image[j_r, j_c] = (
                            image[j_r, j_c]  # remove blank
                            # LM201103 fixed bug renmoved -1 from th
                            # + reference_image[(j_r-1)*nrowbin+j_br+nrowskip-1,(j_c-1)*ncolbinC*ncolbinF+j_bc+ncolskip-1] #row and column value evaluation, -1 to adjust for python indexing
                            + reference_image[
                                (j_r) * nrowbin + j_br + nrowskip,
                                (j_c) * ncolbinC * ncolbinF + j_bc + ncolskip,
                            ]  # row and column value evaluation
                        )

    binned_image = image / gain

    return binned_image, header


def bin_image_using_predict_and_get_true_image(header, reference_image="999"):
    simage_raw_binned, header = bin_image_using_predict(header, reference_image)
    simage_raw_binned, _ = get_true_image(header, simage_raw_binned)
    return simage_raw_binned
