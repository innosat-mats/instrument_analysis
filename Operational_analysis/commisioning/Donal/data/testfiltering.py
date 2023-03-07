#%%
from pyarrow import fs, schema, string
import pyarrow.dataset as ds
import boto3
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from datetime import datetime, timezone
import numpy as np
import matplotlib.pylab as plt
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
import io
from tangentlib import *
#%matplotlib widget
#%%
session = boto3.session.Session(profile_name="mats")
credentials = session.get_credentials()

s3 = fs.S3FileSystem(
    secret_key=credentials.secret_key,
    access_key=credentials.access_key,
    region=session.region_name,
    session_token=credentials.token)
#%%
dataset = ds.dataset(
    "ops-payload-level1a-v0.5",
    filesystem=s3   
)
table = dataset.to_table(
    filter=(
        ds.field('EXPDate') > DT.datetime(2023, 2, 21, 00, 00, tzinfo=timezone.utc)
    ) & ( 
        ds.field('EXPDate') < DT.datetime(2023, 2, 21, 15, 55, tzinfo=timezone.utc)
    ) & (ds.field('TPlat') >-20) & (ds.field('TPlat') < 20)
)
#%%
df = table.to_pandas().reset_index().set_index('TMHeaderTime')
# %%
df
# %%
plt.figure()
plt.plot(df2.TPlon,df2.TPlat,'.')
# %%
df.columns
# %%
startdate=DT.datetime(2023, 2, 21, 00, 00, tzinfo=timezone.utc)
stopdate=DT.datetime(2023, 2, 21, 15, 55, tzinfo=timezone.utc)
#filter=(ds.field('TPlat') >-20) & (ds.field('TPlat') < 20)
filter={"TPlat":[-20,20],"TPlon":[0,90]}
df2=read_MATS_data(startdate,stopdate,version='0.5')
# %%
