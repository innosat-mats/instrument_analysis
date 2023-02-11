#%%
from pyarrow import fs, schema, string
from pyarrow.dataset import FilenamePartitioning
import pyarrow.dataset as ds
import boto3
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pylab as plt
from datetime import datetime,timezone
from scipy.spatial.transform import Rotation as R
from tangentlib import *
from mats_utils.geolocation.coordinates import *
%matplotlib widget
#%%
session = boto3.session.Session(profile_name="mats")
credentials = session.get_credentials()

s3 = fs.S3FileSystem(
    secret_key=credentials.secret_key,
    access_key=credentials.access_key,
    region=session.region_name,
    session_token=credentials.token)

# %%
dataset = ds.dataset(
    "ops-payload-level1b-v0.3/ops-payload-level1a-v0.4",
    filesystem=s3   
)
table = dataset.to_table(
    filter=(
        ds.field('EXPDate') > datetime(2023, 2, 4, 0, 0, tzinfo=timezone.utc)
    ) & ( 
        ds.field('EXPDate') < datetime(2023, 2, 4, 11, 55, tzinfo=timezone.utc)
    ) 
)
# %%
df = table.to_pandas()
# %%
test=xr.Dataset(df)



# %%
plt.figure(figsize=(8,2))
plt.pcolor(np.stack(test.ImageCalibrated[100].values.item()))
plt.colorbar()
plt.clim([0, 20])
plt.show()
#
# %%
plt.figure()
pl
# %%
test.ImageCalibrated[10]
# %%
for var in test.variables:
    print (var,test[var][0])
# %%
