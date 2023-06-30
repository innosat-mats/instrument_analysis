#%%
import numpy as np
import sqlite3 as sqlite
from glob import glob
import os
import pickle
import datetime as DT
#%%
homecat = os.environ.get('HOME')
files= glob(homecat + '/Downloads/hpms/*')
#%%

                     
# %%
def gethpm (date,channelname):
    """
    Function to get the hotpixel map for a given date

    Arguments
    ----------
    date : datetime item specifing the desired date
    channelname : The name of the channel (eg 'IR1') for which the map is required

    Returns
    -------
    mapdate : datetime item giving the date of the map
        if no valid map this will be the same as the date requested 
        
    HPM : array[unit16] or empty array if no valid data
        map of hotpixel counts for the given date 
    """

    db = sqlite.connect(homecat + '/Downloads/hpms/hpms.db')
    cur = db.cursor()
    selectstr= 'select datetime, HPM from hotpixelmaps WHERE  datetime <= "{}"  and channel ==  "{}"  ORDER BY datetime DESC limit 1'
    cur.execute(selectstr.format (date,channelname))
    row=cur.fetchall()
    if row :
        row=row[0]
        mapdate=row[0]
        HPM=pickle.loads(row[1])
    else:
        mapdate= date
        HPM=np.array([])
    return (mapdate,HPM)
# %%
testdate=DT.datetime(2023,2,24)
mapdate,HPM=gethpm(testdate,'IR1')
print (mapdate, HPM)
# %%
