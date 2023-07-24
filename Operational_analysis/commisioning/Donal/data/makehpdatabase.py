#%%
import numpy as np
import sqlite3 as sqlite
from glob import glob
import os
import pickle
import datetime as DT
#%%
homecat = os.environ.get('HOME')
files= glob(homecat + '/Downloads/hpms-2/*')
#%%
db = sqlite.connect(homecat + '/Downloads/hpms/hpms.db')
cur = db.cursor()
cur.execute("create table if not exists hotpixelmaps ('key' INTEGER not null,'datetime' REAL ,'channel' TEXT,'HPM' blob, PRIMARY KEY (key))")
insertstr = "insert or replace into hotpixelmaps values (?, ?, ?, ? )"

# %%
for f in files:
    filename= f.split("/")[-1]
    ch=filename[0:3]
    date=DT.datetime.strptime( filename[4:-4], "%Y%m%d%H%M%S")
    print(ch,date)
    data=np.fromfile(f, dtype=np.uint16)
    hot=data.reshape(187,44)     
    datekey=date.year*100000000+date.month*1000000+date.day*10000+date.hour*100+date.minute
    pic=pickle.dumps(hot, protocol=pickle.HIGHEST_PROTOCOL)                
    cur.execute(insertstr,(datekey,date,ch,pic))
db.commit()
                     
# %%
