#%%
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage, calibrate_CCDitems
import datetime as DT
from selection_tools.itemselect import select_on_time as seltime
from rawdata.time_tools import add_datetime as add_datetime
from geolocation import satellite as satellite
from imagetools.imagetools import shift_image
from plotting.plotCCD import orbit_plot, simple_plot
import pandas as pd

def collapsandplot(imagecub,collapsdim, fix, ax, signallabel='', title=''):
    img_hmean=np.nanmean(imagecube, axis=collapsdim)
    img_mean=np.nanmean(img_hmean, axis=0)
    myindex1=np.arange(0, img_mean.shape[0])
    

    if collapsdim==2:
        myindex=myindex1*CCDitem['NRBIN']#/7.6+54
        ax.plot(img_mean,myindex, label='mean')

        img_std=img_hmean.std(0)
        ax.plot(img_mean+img_std,myindex, '--', label='mean+1std')
        ax.plot(img_mean-img_std,myindex, '--', label='mean+1std')


        ax.plot(img_hmean.max(0),myindex, '.', label='max')
        ax.plot(img_hmean.min(0),myindex, '.', label='min')

        ax.set_xlabel(signallabel)
        ax.set_ylabel('bin nr')
        ax.set_title(title)
        ax.legend()
        plt.tight_layout()
    elif collapsdim==1:
        myindex=myindex1*CCDitem['NCBIN CCDColumns']#/7.6-125
        ax.plot(myindex, img_mean, label='mean')

        img_std=img_hmean.std(0)
        ax.plot(myindex,img_mean+img_std, '--', label='mean+1std' )
        ax.plot(myindex,img_mean-img_std, '--', label='mean-1std')


        ax.plot(myindex,img_hmean.max(0), '.', label='max')
        ax.plot(myindex,img_hmean.min(0), '.', label='min')

        ax.set_ylabel(signallabel)
        ax.set_xlabel('bin nr')
        ax.set_title(title)
        
        ax.legend()
        plt.tight_layout()

    else:
        raise Warning('collapsdim must be 1 or 2')

    return

def create_imagecube(CCDitems, image_specification):
    """    
    Parameters
    ----------
    CCDitems : LIST of CCDitems
    image_specification STRING that can be either 'IMAGE' or 'image_calibrated'

    Returns
    -------
    3d ndarray with all images and time as the last dimension

    """
    imagelist=[]
    if image_specification!='IMAGE' and image_specification!='image_calibrated' and image_specification!='image_common_fov':
        raise KeyError('image specification must be IMAGE or image_alibrated or image_common_fov')
    for CCDitem in CCDitems:
        image=CCDitem[image_specification]
        imagelist.append(image)

    imagecube=np.array(imagelist)

    return imagecube



instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/'
RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_all_old/'
#RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_from25nov/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'

_,df = read_CCDdata(RacOut)
df_all=df.copy()


# %%

df=df_all


date1 = DT.datetime(2022,11,23,10,00,00)
date2 = DT.datetime(2022,11,23,14,00,00)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022


CCDSELECT=6


df = df[(df.CCDSEL == CCDSELECT) & 
(df['NROW'] ==511) &
(df['NCOL'] ==2047) &
(df.TEXPMS ==16000)]


#df[(df.CCDSEL == CCDSELECT)]['TEXPMS'].hist()

#df=df.iloc[0::10, :]


#%%
CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))

# #%%
# CCDitemsdf = pd.DataFrame.from_dict(CCDitems)
# simple_plot(CCDitemsdf, image_path, nstd=2, cmap='inferno', custom_cbar=False,
#     ranges=[0, 1000], format='png')



#%%

calibrate=False
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:
        totbin=CCDitem['NCBIN FPGAColumns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin

if calibrate:
    signallabel='Signal [10^10*ph/cm2/str/nm]'
    image_specification='image_calibrated'
else:
    signallabel='Counts'
    image_specification='IMAGE'

for CCDitem in CCDitems:
    #Shift image, i.e. put image on common field of view
    image_common_fov, error_flags_flipnshift = shift_image(CCDitem, CCDitem[image_specification])
    CCDitem['image_common_fov']=image_common_fov



image_specification='IMAGE'
for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    image=CCDitem[image_specification]
    sp=plot_CCDimage(image, fig, ax, title=CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000))



imagecube=create_imagecube(CCDitems, 'image_common_fov')



#%%

fig, ax = plt.subplots(2,1)

#Horizontal mean
CCDitem=CCDitems[0]
collapsandplot(imagecube,2, fig, ax[0], signallabel=signallabel,
    title=CCDitem['channel']+' EXPT: '+str(CCDitem['TEXPMS'])+ ' NCBIN: '
    +str(CCDitem['NCBIN CCDColumns'])+' NRBIN: '+str(CCDitem['NRBIN']))

#Vertical mean
collapsandplot(imagecube,1,fig, ax[1], signallabel=signallabel,
    title=CCDitem['channel']+' EXPT: '+str(CCDitem['TEXPMS'])+ ' NCBIN: '
    +str(CCDitem['NCBIN CCDColumns'])+' NRBIN: '+str(CCDitem['NRBIN']))

#ax[1].set_xlabel('approx km')
#ax[0].set_ylabel('approx km')
#plt.show()

plt.savefig(image_path+'OperationalAnalysisCropUnbinned_'+CCDitem['channel']+'.png', dpi=700)


# %%
