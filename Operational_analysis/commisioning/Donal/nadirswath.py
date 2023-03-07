# %%
from pyarrow import fs, schema, string
import pyarrow.dataset as ds
import boto3
import numpy as np
import pandas as pd
from datetime import datetime, timezone,timedelta
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from tangentlib import *
import mats_utils.geolocation.coordinates as coordinates
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
%matplotlib widget

# %%
session = boto3.session.Session(profile_name="mats")
credentials = session.get_credentials()

s3 = fs.S3FileSystem(
    secret_key=credentials.secret_key,
    access_key=credentials.access_key,
    region=session.region_name,
    session_token=credentials.token)

# %%
dataset = ds.dataset(
    "ops-platform-level1a-v0.3/ReconstructedData",
    filesystem=s3,
    )

table = dataset.to_table(
    filter=(
        ds.field('time') > pd.to_datetime('2023-2-28T0:0:0z').to_datetime64()
    ) & (
        ds.field('time') < pd.to_datetime('2023-3-01T0:0z').to_datetime64()
    )
)

df = table.to_pandas()
df

# %%
q=np.vstack(df.afsAttitudeState)

# %%
npoints=len(df)
instaxis=np.array([0, 0, -1])
metoOHB  = R.from_matrix([[0,0,-1],[0,-1,0],[-1,0,0]])
nadcam = R.from_euler('XYZ', [0, -(90-225),0], degrees=True).apply([1, 0, 0])
ts=sfapi.load.timescale()

def funheight(s, t, pos, FOV):
    newp = pos + s * FOV
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    return wgs84.subpoint(newp).elevation.m**2


def findsurface(t, pos, FOV):
    res = minimize_scalar(funheight, args=(t, pos, FOV), bracket=(3e5, 8e5))
    return res

#%%
heights=np.zeros(npoints)
satheights=np.zeros(npoints)
points=[]
for i in range (npoints):
    t=ts.from_datetime(df.time[i].replace(tzinfo=sfapi.utc))   
    quat=R.from_quat(np.roll(q[i,:],-1))
    FOVn=quat.apply(metoOHB.apply(nadcam))
    pos=df.afsGnssStateJ2000[i][0:3]
    satheights[i]=wgs84.subpoint(ICRF(Distance(m=pos).au, t=t, center=399)).elevation.m
    res=findsurface(t,pos,FOVn)
    heights[i]=res.x
    newp = pos + res.x * FOVn
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    points.append(newp)


# %%
plt.close('all')
plt.figure(figsize=[7,5])
plt.plot(RAs,Decs,'+')
plt.plot(RAs[0:15],Decs[0:15],'r+')
plt.plot(*xyz2radec(np.expand_dims(meanvec,axis=1),deg=True),'bo')
plt.xlabel('Ra (deg)')
plt.ylabel('Dec (deg)')

# %%

plt.figure(figsize=(6,3))


plt.plot(df.time[:npoints],RAs)
ax1=plt.gca()
ax2 = ax1.twinx()
ax2.plot(df.time[:npoints],Decs,'r.')
plt.gcf().autofmt_xdate()
ax1.set_ylabel('Ra')
ax2.set_ylabel('Dec', color='r')

# %%
mask=(df.time>datetime(2022,12,30,6,15,23)) & (df.time<datetime(2022,12,30,6,15,23)+timedelta(milliseconds=15000))
np.vstack(df[mask].afsAttitudeState).mean(axis=0)

# %%
coordinates.meanquaternion(datetime(2022,12,30,6,15,23),timedelta(milliseconds=15000))

# %%
plt.figure()
plt.hist(angles,list(np.linspace(0,0.35,18)))
plt.ylabel('Number of points')
plt.xlabel('Arc minutes')

# %%
pixelmap=np.zeros((21,21))
dispersion=6.06 / 2047*60   # arcmin/pix

# %%
dispersion

# %%
np.linspace(0,0.35,36)

# %%

dataset = ds.dataset(
    "ops-platform-level1a-v0.3/ReconstructedData",
    filesystem=s3,
    #ignore_prefixes=['PreciseAttitude', 'HK', 'TM', 'sco']
    #ignore_prefixes=['PreciseOrbitEstimation', 'HK', 'TM', 'sco']

)

table = dataset.to_table(
    filter=(
        ds.field('time') > pd.to_datetime('2023-2-13T12:00z').to_datetime64()
    ) & (
        ds.field('time') < pd.to_datetime('2023-2-15T12:00z').to_datetime64()
    )
)



# %%
df = table.to_pandas()

# %%
plt.figure()

plt.plot(df.time,df.afsTangentH_wgs84)
plt.gcf().autofmt_xdate()

# %%
df

# %%



