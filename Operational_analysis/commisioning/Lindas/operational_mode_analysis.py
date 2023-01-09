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


def collapsandplot(imagecub,collapsdim, fix, ax, signallabel='', title=''):
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

def create_imagecube(CCDitems, image_specfication='IMAGE'):
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
RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_all/'
#RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_from25nov/'
calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'

_,df = read_CCDdata(RacOut)
df_all=df.copy()


# %%

df=df_all

#hej=df_all[(df.NRBIN==200)]

date1 = DT.datetime(2022,12,1,12,00,00)
date2 = DT.datetime(2022,12,30,12,00,00)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022



#CCD limb settings:
CCDSELECT=6
df[(df.CCDSEL == CCDSELECT)]['NCBIN CCDColumns'].hist(bins=100)
#%%
if (CCDSELECT==2 or CCDSELECT==3): #IR4 or IR3
    NRBINop=6
    NCBIN_CCDop=200
    NCBIN_FPGAop=1
    NROWop=85
    NCOLop=9
    TEXPMS=1500
elif (CCDSELECT==1 or CCDSELECT==4): #IR1 or IR2
    NRBINop=2
    NCBIN_CCDop=40
    NCBIN_FPGAop=1
    NROWop=255
    NCOLop=50
    TEXPMS=5000
elif (CCDSELECT==5 or CCDSELECT==6): #UV1 or UV2
    NRBINop=2
    NCBIN_CCDop=40
    NCBIN_FPGAop=1
    NROWop=255
    NCOLop=50
    TEXPMS=5000
elif (CCDSELECT==7): # NADIR
    NRBINop=36
    NCBIN_CCDop=36
    NCBIN_FPGAop=1
    NROWop=14
    NCOLop=55
    TEXPMS=2000
else: 
    raise Warning('undefined CCD')


df = df[(df.CCDSEL == CCDSELECT) & 
(df.NRBIN ==NRBINop)& 
(df['NCBIN CCDColumns'] ==NCBIN_CCDop) &
(df['NCBIN FPGAColumns'] ==NCBIN_FPGAop) &
(df['NROW'] ==NROWop) &
(df['NCOL'] ==NCOLop) &
(df.TEXPMS ==TEXPMS)]
#%%

#df=df.iloc[0::10, :]

CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))



calibrate=False
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:
        totbin=CCDitem['NCBINFPGA Columns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin

if calibrate:
    signallabel='Signal [10^10*ph/cm2/str/nm]'
    image_specification='image_calibrated'
else:
    signallabel='Counts'
    image_specification='IMAGE'

CCDitems_night=[]
CCDitems_day=[]
for CCDitem in CCDitems:
    # Select day or night 
    (satlat, satlon, satLT, nadir_sza, nadir_mza,
        TPlat, TPlon,TPLT, TPsza, TPssa) = satellite.get_position(CCDitem['EXP Date'])
    #if TPsza >90: #Nighttime
    if TPlat <-60: #Nighttime    
        CCDitems_night.append(CCDitem)
    else:
        CCDitems_day.append(CCDitem)

#%%
dayornight='night'
if dayornight=='day': 
    imagecube=create_imagecube(CCDitems_day, image_specification)
else:
    imagecube=create_imagecube(CCDitems_night, image_specification)

#%%

fig, ax = plt.subplots(2,1)
if calibrate:
    signallabel='Signal [10^10*ph/cm2/str/nm]'
else:
    signallabel='Counts'

CCDitem=CCDitems[0]
#Horizontal mean
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

plt.savefig(image_path+'OperationalAnalysis_saturation_'+CCDitem['channel']+'.png', dpi=700)

#%%

values=imagecube.mean(0).mean(1)
airglowstartheitht=100
a=values[airglowstartheitht:250]
maxind=np.where(a==np.max(a))
ind=maxind[0][0]+airglowstartheitht
imagecube.mean(0).mean(1)[ind]


# Creating histogram
NLCmaxUV2=215
fig2, ax2 = plt.subplots(figsize =(10, 7))
ax2.hist(imagecube[:,NLCmaxUV2,10:40].mean(1), bins = 100)
plt.title(dayornight+' at NLC height')
# %%
