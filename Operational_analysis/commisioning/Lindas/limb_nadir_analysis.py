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
#from plotting.plotting_functions import collapsandplot, create_imagecube
from rawdata.time_tools import add_datetime as add_datetime
from geolocation import satellite as satellite

def collapsandplot(imagecub,collapsdim, ax, signallabel='', title=''):
    """    
    Parameters
    ----------
    imagecub : 3d numpy array made by several CCDimages, created by  create_imagecube(CCDitems, calibrated=False)
    calibrated : BOOLEAN 

    Returns
    -------
    3d ndarray with all images and time as the last dimension

    """

    img_hmean=imagecube.mean(collapsdim)
    img_mean=img_hmean.mean(0)
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

def create_imagecube(CCDitems, calibrated=False):
    """    
    Parameters
    ----------
    CCDitems : LIST of CCDitems
    calibrated : BOOLEAN 

    Returns
    -------
    3d ndarray with all images and time as the last dimension

    """
    imagelist=[]
    for CCDitem in CCDitems:
        image=CCDitem[image_specification]
        imagelist.append(image)

    imagecube=np.array(imagelist)

    return imagecube



instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/'
RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_LimbPointingInNadir/'
#RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_from25nov/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'

_,df = read_CCDdata(RacOut)
df_all=df.copy()

df[(df.CCDSEL == 3)]['TEXPMS'].hist()
# %%

df=df_all

#hej=df_all[(df.NRBIN==200)]

# date1 = DT.datetime(2022,11,23,12,00,00)
# date2 = DT.datetime(2022,11,30,12,00,00)

# dfdatetime=add_datetime(df)
# df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022



#CCD limb settings:
CCDSELECT=1


df = df[(df.CCDSEL == CCDSELECT) ]


#%%

#df=df.iloc[0::100, :]
#%%
CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))

#%%

calibrate=False
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:
        totbin=CCDitem['NCBINFPGA Columns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin



#%%

for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    image=CCDitem[image_specification]
    sp=plot_CCDimage(image, fig, ax, title=CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000))


#%%
#plt.savefig(image_path+'OperationalAnalysis2'+CCDitem['channel']+'.png', dpi=1600)

# %%
