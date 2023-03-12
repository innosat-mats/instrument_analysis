#%%
import numpy as np
import os
from datetime import datetime, timezone,timedelta
import matplotlib.pylab as plt
import matplotlib as mpl
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve
from tangentlib import *
import mats_utils.geolocation.coordinates as coordinates
from skyfield import api as sfapi
import skyfield.sgp4lib as sgp4lib
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
#from getmatstle import *
#%%
ts=sfapi.load.timescale()
planets = sfapi.load('de421.bsp')
earth=planets['Earth']
sun=planets['Sun']
def makeplot(number,starttime):
    deltat=timedelta(seconds=0.5*60)
    #starttime=datetime(2023,3,1,0,0,tzinfo=timezone.utc)
    times=[starttime + n * deltat for n in range(200)]
    #tle = get_tle_dateDB(date,8)
    tle=['1 26702U 01007A   23037.91034836  .00007274  00000+0  42696-3 0  9993', '2 26702  97.5120  53.5306 0009881 191.0534 169.0481 15.11989832201676']
    sfmats = sgp4lib.EarthSatellite(tle[0],tle[1])
    period= 2*np.pi/sfmats.model.nm
    t=ts.from_datetimes(times)
    g=sfmats.at(t)
    altaz=(earth+wgs84.subpoint(g)).at(t).observe(sun).apparent().altaz()
    Nadirsza=[(90-ang) for ang in altaz[0].degrees]
    sublats=g.subpoint().latitude.degrees
    sublons=g.subpoint().longitude.degrees
    ECI_pos=g.position.m
    ECI_vel=g.velocity.m_per_s
    vunit=np.array(ECI_vel)/norm(ECI_vel,axis=0)
    mrunit=-np.array(ECI_pos)/norm(ECI_pos,axis=0)
    yunit=np.cross(mrunit[:,1],vunit[:,1]) #constant arround the orbit
    #print (sublat,sublon)
    ascnod=-np.cross([0,0,1],yunit)
    ascnod/=norm(ascnod)
    arg=np.array([np.rad2deg(np.arccos(np.dot(-mrunit[:,i],ascnod.T)))*np.sign(-mrunit[2,i]) for i in range(vunit.shape[1])] ).T
    yaw =0
    yawoffset=0
    pitch=-20.86 #deg
    print(pitch)
    #pitch=findpitch(87500,t[1], ECI_pos[:,1], np.deg2rad(yaw)+yawoffset, rotmatrix)
    yaw=-3.3*np.cos(np.deg2rad(arg-pitch+20)) #deg
    yaw=yaw#*0
    #print (np.rad2deg(pitch))
    FOV=[]
    for i in range(vunit.shape[1]):
        FOVi = R.from_euler('XYZ', [0, pitch, yaw[i]], degrees=True).apply([1, 0, 0])
        rotmatrix=np.array([vunit[:,i],yunit,mrunit[:,i]]).T 
        #FOV=rotate(np.array([1,0,0]),np.deg2rad(yaw)+yawoffset,pitch,0,deg=False)
        FOV.append(np.matmul(rotmatrix,FOVi))
    #res = findtangent(t[1],ECI_pos[:,1],FOV)
    #s=res.x
    s=24535080 # make a fixed value for this
    print('s = ', s)
    FOV=np.array(FOV).T
    newp = ECI_pos + s * FOV
    #print (newp)
    #    pos_s=np.matmul(itrs.rotation_at(t),newp)
    newp=ICRF(Distance(m=newp).au,t=t,center=399)
    tplats=(wgs84.subpoint(newp).latitude.degrees)
    tplons=(wgs84.subpoint(newp).longitude.degrees)
    sundir=(earth+wgs84.subpoint(newp)).at(t).observe(sun).apparent()
    altaz=sundir.altaz()
    TPsza=[(90-ang) for ang in altaz[0].degrees]
    TPssa = (np.rad2deg(np.arccos(np.dot(FOV.T,(sundir.position.m/norm(sundir.position.m,axis=0))))))[:,0]

    #FOV=np.matmul(rotmatrix,FOV)
    # v_dash=FOV
    # y_dash=-np.cross(v_dash,mrunit)
    # y_dash/=norm(y_dash)
    # [FOV_ra,FOV_dec]=xyz2radec(v_dash,deg=True)
    # r_dash=np.cross(v_dash,y_dash)
    # r_dash/=norm(r_dash)
    # invrotmatrix=np.linalg.inv(np.array([v_dash,y_dash,r_dash]).T) 

    def add_arrows(x,y,clinecolour):
        for i in range(len(y)-1):
            if y[i]*y[i+1] < 0 : #sign change        
                plt.arrow(x[i], y[i], x[i+1]-x[i], y[i+1]-y[i], 
                        shape='full', lw=0, length_includes_head=True, head_width=2,color=clinecolour)
    if not os.path.isdir('/tmp/nadirpics'):
        os.mkdir('/tmp/nadirpics')
    plt.figure(1)
    plt.clf()
    plt.plot(Nadirsza,sublats)
    plt.xlim([50,130])
    plt.ylim([-90,90])    
    plt.ylabel('Latitude')
    plt.xlabel('Solar zenith angle')
    plt.title('Nadir Camera SZA  DOY {}'.format(times[0].timetuple().tm_yday) )
    clinecolour=plt.gca().get_lines()[0].get_color()
    add_arrows(Nadirsza,sublats,clinecolour)
    filename='/tmp/nadirpics/Nadir{:04d}.png'.format(number)
    plt.savefig(filename)

    if not os.path.isdir('/tmp/limbpics'):
        os.mkdir('/tmp/limbpics')
    plt.figure(2)
    plt.clf()
    plt.plot(TPsza,tplats,TPssa,tplats)
    plt.xlim([50,130])
    plt.ylim([-90,90])    
    plt.ylabel('Latitude')
    plt.xlabel('Solar zenith angle')
    plt.title('Limb Tangent point SZA/SSA  DOY {}'.format(times[0].timetuple().tm_yday) )
    plt.legend(['SZA','SSA'])
    clinecolour=plt.gca().get_lines()[0].get_color()
    add_arrows(TPsza,tplats,clinecolour)
    clinecolour=plt.gca().get_lines()[1].get_color()
    add_arrows(TPsza,tplats,clinecolour)
    filename='/tmp/nadirpics/Limb{:04d}.png'.format(number)
    plt.savefig(filename)
    return yaw
        
# %%
steps=timedelta(hours=24)
starttime=datetime(2023,1,1,0,0,tzinfo=timezone.utc)
dates=[starttime + n* steps for n in range(365)]
for i in range(365):
    yaw = makeplot(i,dates[i])
# %%
