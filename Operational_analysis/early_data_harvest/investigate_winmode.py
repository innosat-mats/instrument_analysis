# 22-10-11 - Early Harvest Analysis - Björn
# Investigating how the automatic selection of WDW picks mode and if OK.
# This concerns pass 72

import xarray as xr
from mats_l1_processing import read_in_functions as rif
import pandas as pd
import matplotlib.pyplot as plt

directory='/home/waves/projects/MATS/data/CommisioningData/EarlyHarvestPlots/Pass72/pass72/'
items = rif.read_CCDitems(directory)
keys = items[0].keys()

# dictonary to pandas
items = pd.DataFrame.from_dict(items)

# Then isolate 'WDW Mode' = 'Automatic'
items_manual = items[items['WDW Mode'] == 'Manual']
items = items[items['WDW Mode'] == 'Automatic']

# Find 'WDW InputDataWindow' for each CCD
#items['WDW InputDataWindow']
WDWs_manual = ['15..0']
WDWs = ['11..0','12..1','13..2','14..3','15..4']

for CCDno in range(1,8):

    fig, ax = plt.subplots(2,3, figsize=(12,8))
    ax = ax.ravel()
    
    CCDs = items[items['CCDSEL'] == CCDno]
    CCDs = CCDs[CCDs['TEXPMS'] == 6000]
    CCDs_man = items_manual[items_manual['CCDSEL'] == CCDno]
    CCDs_man = CCDs_man[CCDs_man['TEXPMS'] == 6000]
    #print(CCDs['WDW InputDataWindow'])
    #ax[i].title(f'CCDSEL {CCDno}')

    i = 0
    for WDWno in WDWs:

        WDW_ccd = CCDs[CCDs['WDW InputDataWindow'] == WDWno]

        for index,row in WDW_ccd.iterrows():
            image = row['IMAGE']
            ax[i].hist(image.flatten(), bins=200, label=f'{row["TEXPMS"]}', alpha=0.6, density=True)
        
        ax[i].set_title(f"{WDWno} (x{len(WDW_ccd['WDW InputDataWindow'])})")
        ax[i].legend(title='TEXPMS')
        i = i + 1
    
    WDW_ccd_man = CCDs_man[CCDs_man['WDW InputDataWindow'] == '15..0']

    for index,row in WDW_ccd_man.iterrows():
            image = row['IMAGE']
            ax[-1].hist(image.flatten(), bins=200, label=f'{row["TEXPMS"]}', alpha=0.6, density=True)
            ax[-1].legend(title='TEXPMS')

    fig.suptitle(f'CCDSEL : {CCDno}')
    plt.tight_layout()

    plt.savefig(f'/home/waves/projects/MATS/scripts/plots/all_wdw/histo_CCD{CCDno}_all_WDW.png', format='png')

#plt.figure()
#CCDs['IMAGE'].plot.hist(stacked=True, bins=20)
# Look at values - Histogram animation / Max / Min / STD / Saturation?

# We are going to generate subplots; 
# For EACH CCD - CCD 1 to 7.







