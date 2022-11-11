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