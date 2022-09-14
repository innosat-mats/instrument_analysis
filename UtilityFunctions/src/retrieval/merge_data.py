import os
import numpy as np
import xarray as xr
from halo import Halo

# settings
nx = 500
overlap = 100
end = int(nx - overlap/2)

# folder containing patches to merge
dir = '/home/waves/projects/hiamcm-juwave/data/fsps/test_orbits_atm/MATS_avg_orbit_data/patches/'
outdir = '/home/waves/projects/hiamcm-juwave/data/fsps/test_orbits_atm/MATS_avg_orbit_data/'

# spinners 
list_spin = Halo(text='loading list of files...', spinner='dots')
save_spin = Halo(text='saving ...', spinner='dots')

def rename_dim(data, dim1 = 'x_tp', dim2 = 'x'):
    # function to change sat_var dim name
    # to avoid concat errors

    satvars = ("SAT_ALT", "SAT_LAT",
               "SAT_LON", "SAT_TIME",
               "TANGENT_POINT_ALT",
               "TANGENT_POINT_LAT",
               "TANGENT_POINT_LON")

    for var in satvars:
        if var in list(data.keys()):
            data[var] = data[var].rename({dim1:dim2})
    return data

list_spin.start()
files =  os.listdir(dir)
files.sort()
list_spin.succeed()

print(f'-------- merging {len(files)} orbit patches ----------')
print(f'----------- approx. {len(files)*400*20 / 40000} orbits -------------')

data = xr.open_dataset(f'{dir}/{files[0]}')
data = data.isel(x=slice(0,end), x_tp=slice(0,end))
data = rename_dim(data)

for i in range(1,len(files)):

    merge_spin = Halo(text=f'merging patch {i+1} ...', spinner='dots')
    merge_spin.start()
    temp = xr.open_dataset(f'{dir}/{files[i]}')
    temp = temp.isel(x=slice(int(overlap/2),end), x_tp=slice(int(overlap/2),end))
    temp = rename_dim(temp)
    data = xr.concat([data, temp], dim='x')
    merge_spin.succeed()

data = rename_dim(data, dim1='x', dim2='x_tp')

# correct lost time dim (from applying kernels)
list_time_dims = ('DENS', 'N2_SMOOTH', 'PRES', 'TEMP', 'TEMP_BACKGROUND','TEMP_RESIDUAL', 'U',
                'U_BACKGROUND', 'U_RESIDUAL', 'V', 'V_BACKGROUND', 'V_RESIDUAL',
                'W', 'W_BACKGROUND', 'W_RESIDUAL')

for dim in list_time_dims:
    data[dim] = data[dim].expand_dims('time')

# transpose to match (t, z, y, x, xtp)
data = data.transpose('time','height','y','x', 'x_tp')

# generate name and save
print(f'-------- all {len(files)} orbit patches merged ----------')
save_spin.start()
outname = files[0].split('_nx')[0]+f'_{i+1}p.nc'
outpath = outdir + outname 
data.load().to_netcdf(path=outpath, mode="w", format="NETCDF4_CLASSIC")
save_spin.succeed()