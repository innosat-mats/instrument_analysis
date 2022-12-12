import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import cartopy.crs as ccrs
import pandas as pd
from geolocation import satellite as satellite
from cartopy.feature.nightshade import Nightshade

def check_type(CCDitems):
    """Check format of CCDitems
    Exit program if type is not DataFrame

    Parameters
    ----------
    CCDitems : any
        CCDitems
    """

    if isinstance(CCDitems, pd.core.frame.DataFrame) is False:
        sys.exit("CCDitems need to be converted to DataFrame!")

    return


def simple_plot(CCDitems, outdir, nstd=2, plot_calibrated=False, cmap='inferno', custom_cbar=False,
                ranges=[0, 1000], format='png'):
    """Generates plots from CCDitems with basic orbit parameters included.
    Images will be sorted in folders based on CCDSEL in directory specified.

    Parameters
    ----------
    CCDitems : DataFrame
        CCDitems for plotting
    outdir : str
        Out directory
    nstd : int, optional
        number of standard deviations, by default 2
    plot_calibrated : bool, optional
        raw or calibrated image'
    cmap : str, optional
       colormap for plot, by default 'inferno'
    custom_cbar : bool, optional
        if custom cbar set True, by default False
    ranges : list, optional
        limits for custom cbar, by default [0,1000]
    format : str, optional
        format for files, by default 'png'
    """

    check_type(CCDitems)

    if plot_calibrated:
        image_str = 'image_calibrated'
    else:
        image_str = 'IMAGE'

    fig = plt.figure(figsize=(12, 3))

    for CCDno in range(0, 8):
        CCDs = CCDitems[CCDitems['CCDSEL'] == CCDno]

        if len(CCDs) > 0:
            outpath = f"{outdir}CCDSEL{str(CCDno)}"
            if not os.path.exists(outpath):
                os.makedirs(outpath)

        for index, CCD in CCDs.iterrows():

            # save parameters for plot
            channel = CCD['channel']
            image = CCD[image_str]
            [col, row] = image.shape
            texpms = CCD['TEXPMS']
            exp_date = CCD['EXP Date'].strftime("%Y-%m-%dT%H:%M:%S:%f")

            # filename
            outname = f"{CCD['Image File Name'][:-4]}_{index}"

            # calc std and mean
            mean = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].mean()
            std = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].std()


            # orbital parameters
            (satlat, satlon, satLT,
             nadir_sza, nadir_mza,
             TPlat, TPlon,
             TPLT, TPsza, TPssa) = satellite.get_position(CCD['EXP Date'])

            # plot CCD image
            if custom_cbar:
                plt.pcolormesh(image, cmap=cmap,
                               vmax=ranges[1], vmin=ranges[0])
            else:
                plt.pcolormesh(image, cmap='magma',
                               vmax=mean+nstd*std, vmin=mean-nstd*std)
            plt.colorbar(label='counts')

            # print out additional information
            plt.figtext(0.1, 0.8, f'tpSZA: {TPsza:.6}',
                        fontsize=10, color='white')
            plt.figtext(0.5, 0.8, (f'satlat, satlon: ({satlat:.6}' +
                                   f', {satlon:.6})'),
                        fontsize=10, color='white')
            plt.figtext(0.25, 0.8, f'TPlat, TPlon: ({TPlat:.6}, {TPlon:.6})',
                        fontsize=10, color='white')

            plt.title(f'ch: {channel}; time: {exp_date}; TEXPMS: {texpms}')
            plt.tight_layout()

            # save figure
            plt.savefig(f'{outpath}/{outname}.{format}', format=format)
            plt.close()


def orbit_plot(CCDitems, outdir, nstd=2, plot_calibrated=False, cmap='inferno', custom_cbar=False,
               ranges=[0, 1000], format='png'):
    """
       Generates plots from CCD items: image, histogram and map.
       If simple_plot is True, only CCD image is plotted.
       Figures will be saved in subfolders of outdir by CCDSEL.

    Parameters
    ----------
    CCDitems : DataFrame
        CCDitems to be plotted.
    outdir : str
        path where images will be saved
    nstd : int, optional
        Number of standard deviations for cbar and histogram, by default 2
    plot_calibrated : bool, optional
        by default False. plot raw or calibrated image'
    cmap : str, optional
        Colourmap for image, by default 'inferno'
    custom_cbar : bool, optional
        Custom cbar, by default False
    ranges : tuple, optional
        If custom_cbar == True, specify cbar limits, by default (0,1000)
    format : str
        file format for output img
    """

    check_type(CCDitems)

    if plot_calibrated:
        image_str = 'image_calibrated'
    else:
        image_str = 'IMAGE'

    for CCDno in range(0, 8):

        CCDs = CCDitems[CCDitems['CCDSEL'] == CCDno]

        if len(CCDs) > 0:
            outpath = f"{outdir}CCDSEL{str(CCDno)}"

            if not os.path.exists(outpath):
                os.makedirs(outpath)

            for index, CCD in CCDs.iterrows():

                # save parameters for plot
                channel = CCD['channel']
                image = CCD[image_str]
                [col, row] = image.shape
                texpms = CCD['TEXPMS']
                exp_date = CCD['EXP Date'].strftime("%Y-%m-%dT%H:%M:%S:%f")

                # filename
                outname = f"{CCD['Image File Name'][:-4]}_{index}"

                # calc std and mean
                mean = image[int(col/2-col*4/10):int(col/2+col*4/10),
                             int(row/2-row*4/10):int(row/2+row*4/10)].mean()
                std = image[int(col/2-col*4/10):int(col/2+col*4/10),
                            int(row/2-row*4/10):int(row/2+row*4/10)].std()

                # orbital parameters
                (satlat, satlon,
                 satLT, nadir_sza,
                 nadir_mza, TPlat,
                 TPlon, TPLT, TPsza,
                 TPssa) = satellite.get_position(CCD['EXP Date'])

                fig = plt.figure(figsize=(10, 7))

                # generate subplot grid
                ax = plt.subplot2grid((2, 2), (1, 0),
                                      colspan=1, rowspan=1,
                                      projection=ccrs.PlateCarree(),
                                      fig=fig)
                ax1 = plt.subplot2grid((2, 2), (0, 0),
                                       rowspan=1, colspan=2, fig=fig)
                ax2 = plt.subplot2grid((2, 2), (1, 1), rowspan=1,
                                       colspan=1, fig=fig)

                # map settings
                gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                                  linewidth=0.9, color='black',
                                  alpha=0.5, linestyle='-')
                gl.xlabels_top = True
                gl.ylabels_left = True
                gl.ylabels_right = False
                gl.xlines = True
                ax.set_xlabel('longitude [deg]')
                ax.set_ylabel('latitude [deg]')
                ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
                ax.add_feature(Nightshade(CCD['EXP Date'], alpha=0.2))
                ax.coastlines()

                # plot CCD image
                if custom_cbar:
                    img = ax1.pcolormesh(image, cmap=cmap,
                                         vmax=ranges[1], vmin=ranges[0])
                else:
                    img = ax1.pcolormesh(image, cmap='magma',
                                         vmax=mean+nstd*std,
                                         vmin=mean-nstd*std)
                ax1.set_title(f'ch: {channel}; time: '
                              + f'{exp_date}; TEXPMS: {texpms}')
                fig.colorbar(img, ax=ax1)

                # plot sat position and tangent point
                ax.scatter(satlon, satlat, s=10,
                           color='red', label='satellite pos.')
                ax.scatter(TPlon, TPlat, s=10,
                           color='green', label='TP pos.')
                ax.legend(ncol=2, fontsize=7, loc='lower right')

                # print out additional information
                plt.figtext(0.15, 0.03, f'nadirSZA: {nadir_sza:.6}',
                            fontsize=10)
                plt.figtext(0.15, 0.06, f'nadirMZA: {nadir_mza:.6}',
                            fontsize=10)
                plt.figtext(0.35, 0.03, f'tpSZA: {TPsza:.6}', fontsize=10)
                plt.figtext(0.35, 0.06, f'tpSSA: {TPssa:.6}', fontsize=10)

                # plot histogram
                nbins = int(1 + np.ceil(np.log2(len(image.flatten()))))
                ax2.hist(image.flatten(), bins=nbins, alpha=0.6,
                         density=True, range=[mean-nstd*std, mean+nstd*std])
                ax2.set_xlabel('counts')
                ax2.axvline(x=mean, label='mean',
                            linestyle='--', linewidth=1.5)
                ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
                ax2.legend(loc='upper right')
                ax2.grid()
                plt.savefig(f'{outpath}/{outname}.{format}', format=format)
                plt.close()

    return
