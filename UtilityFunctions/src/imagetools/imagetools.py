import numpy as np

def bin_image(image,nrbin,ncolbin):
    """Function which bins an image

    Arguments:
        image (np.array): image to be binned
        nrbin (int): number of rows to be binned
        ncolbin (int): number of columns to be binned

    Returns:
        (np.array): binned image

    """
    nrow,ncol = image.shape()
    nrow_binned = np.floor(nrow/nrbin)
    ncol_binned = np.floor(ncol/ncolbin)

    colbin = np.zeros([nrow,ncol_binned])
          
    for j in range(0,ncol_binned):
        colbin[:,j] = image[:,j*ncolbin:(j+1)*ncolbin].sum(axis=1)

    # declare zero array for row binning 
    binned = np.zeros([nrow_binned,ncol_binned])
    
    for j in range(0,nrow_binned):
        binned[j,:] = colbin[j*nrbin:(j+1)*nrbin,:].sum(axis=0)

    return binned

def shift_image(CCDitem, image=None):
    """ 
    Shift the images to account for the misalignment. 
    Or rather put the image on a common field of view with all other channels.
    Args:
        CCDitem
        optional image 

    Returns: 
        
        image that has been shifted
        error_flag

    """
    from mats_l1_processing.grid_image import get_shift

    x_pos,y_pos = get_shift(CCDitem)

    x_maximum=75 #156 This is the maximum shift, silly that it is hardcoded
    y_maximum=192
    x_rel=x_maximum-x_pos
    y_rel=y_maximum-y_pos
    
    
    image_common_fov = np.empty((720,2300))
    error_flag= np.ones(image_common_fov.shape, dtype=np.uint16)
    image_common_fov[:] = np.nan
    image_common_fov[y_rel:y_rel+image.shape[0], x_rel:x_rel+image.shape[1]]=image
    error_flag[y_rel:y_rel+image.shape[0], x_rel:x_rel+image.shape[1]] = 0
    
    return image_common_fov, error_flag    