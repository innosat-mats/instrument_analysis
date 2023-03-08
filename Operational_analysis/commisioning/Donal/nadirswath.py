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
        ds.field('time') > pd.to_datetime('2023-1-06T0:0:0z').to_datetime64()
    ) & (
        ds.field('time') < pd.to_datetime('2023-2-07T0:0z').to_datetime64()
    )
)

df = table.to_pandas()
df=df.iloc[range(0,len(df),10)]

# %%
q=np.vstack(df.afsAttitudeState)

# %%
npoints=len(df)
instaxis=np.array([0, 0, -1])
metoOHB  = R.from_matrix([[0,0,-1],[0,-1,0],[-1,0,0]])
nadcam = R.from_euler('XYZ', [0, -(90-23),0], degrees=True).apply([1, 0, 0])
nadcamr= R.from_euler('XYZ', [12.2, -(90-23),0], degrees=True).apply([1, 0, 0])
nadcaml= R.from_euler('XYZ', [-12.2, -(90-23),0], degrees=True).apply([1, 0, 0])
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
cpoints=[]
lpoints=[]
rpoints=[]
for i in range (npoints):
    t=ts.from_datetime(df.time.iloc[i].replace(tzinfo=sfapi.utc))   
    quat=R.from_quat(np.roll(q[i,:],-1))
    FOVn=quat.apply(metoOHB.apply(nadcam))
    FOVl=quat.apply(metoOHB.apply(nadcaml))
    FOVr=quat.apply(metoOHB.apply(nadcamr))
    pos=df.afsGnssStateJ2000.iloc[i][0:3]
    satheights[i]=wgs84.subpoint(ICRF(Distance(m=pos).au, t=t, center=399)).elevation.m
    res=findsurface(t,pos,FOVn)
    heights[i]=res.x
    newp = pos + res.x * FOVn
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    cpoints.append(newp)
    res=findsurface(t,pos,FOVl)
    heights[i]=res.x
    newp = pos + res.x * FOVl
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    lpoints.append(newp)
    res=findsurface(t,pos,FOVr)
    heights[i]=res.x
    newp = pos + res.x * FOVr
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    rpoints.append(newp)




# %%
clats=[wgs84.subpoint(t).latitude.degrees for t in cpoints]
rlats=[wgs84.subpoint(t).latitude.degrees for t in rpoints]
llats=[wgs84.subpoint(t).latitude.degrees for t in lpoints]
clons=[wgs84.subpoint(t).longitude.degrees for t in cpoints]
rlons=[wgs84.subpoint(t).longitude.degrees for t in rpoints]
llons=[wgs84.subpoint(t).longitude.degrees for t in lpoints]
# %%
planets = sfapi.load('de421.bsp')
earth=planets['Earth']
sun=planets['Sun']
npoints=600
szal=np.zeros(npoints)
szac=np.zeros(npoints)
szar=np.zeros(npoints)
for i in range(npoints):
    t=ts.from_datetime(df.time.iloc[i].replace(tzinfo=sfapi.utc))  
    szal[i]=90-((earth+wgs84.subpoint(lpoints[i])).at(t).observe(sun).apparent().altaz())[0].degrees
    szac[i]=90-((earth+wgs84.subpoint(cpoints[i])).at(t).observe(sun).apparent().altaz())[0].degrees    
    szar[i]=90-((earth+wgs84.subpoint(rpoints[i])).at(t).observe(sun).apparent().altaz())[0].degrees

# %%
plt.figure()
plt.scatter([rlons[0:npoints],clons[0:npoints],llons[0:npoints]],[rlats[0:npoints],clats[0:npoints],llats[0:npoints]],c=[szar,szac,szal],s=2)
plt.colorbar()
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.title('2023-01-06')
# %%
