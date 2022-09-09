# Script for applying averaging kernels to HIAMCM orbit data; 
# %%
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from retrieval.averaging_kernels import apply_3d_kernel
import time
import cProfile
import io
import pstats

# %% OPTIONS

# file information
main_path = "/home/waves/projects/hiamcm-juwave/data/fsps/test_orbits_atm/MATS_orbit_data/no_mid_track/"
orbit_file = "mats_orbit_12t_D5844_97.6_ls0_0600UTZ01JAN2016_spline.nc"
out_path = "/home/waves/projects/hiamcm-juwave/data/fsps/test_orbits_atm/MATS_avg_orbit_data/"

# parallel processing 
parallel = True

# run profiler (requires that parallel is False)
profiler = False

# save as npy instead of xarray
save_npy = False

# fwhm (km)
fwhm_x, fwhm_y, fwhm_z = 80, 5, 1

# along-track points
nx = 320

for nx in [2000]:
    # %% READ DATA
    print('initiating... \n')
    data = xr.open_dataset(main_path + orbit_file)
    data = data.isel(time=0, x=slice(0, nx), x_tp=slice(0, nx))
    TEMP = data.TEMP

    # %% APPLY TO XN CROSS SECTIONS
    TEMP = TEMP.T.values
    x = data.x.values
    y = data.y.values
    z = data.height.values

    if profiler:
        # start profiling
        pr = cProfile.Profile()
        pr.enable()

    runtime = time.time()

    # apply 3D averaging kernels
    print('applying kernels... \n')
    averaged_data = apply_3d_kernel(TEMP, x, y, z, [fwhm_x, fwhm_y, fwhm_z], only_kernel=False, pp=parallel)

    elapsed = time.time() - runtime

    if profiler:
        # end profiling
        pr.disable()
        pr.print_stats(sort='time')

        # save profiler data
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
        ps.print_stats()
        with open("/home/waves/projects/instrument_analysis/UtilityFunctions/src/retrieval/logs/profile.txt", 'w+') as f:
            f.write(s.getvalue())

    # save nx - time data
    with open("/home/waves/projects/instrument_analysis/UtilityFunctions/src/retrieval/logs/time.txt", "ab") as f:
        np.savetxt(f, np.c_[nx, elapsed])

    # convert from list to array
    averaged_data = np.asarray(averaged_data)

    # save data
    print('saving... \n')
    if save_npy:
        np.save(out_path + orbit_file[:-3] + f'_avg_x{fwhm_x}y{fwhm_y}z{fwhm_z}_nx{nx}',
                averaged_data)

    else:
        data.TEMP.values = averaged_data.T
        data.load().to_netcdf(path=(out_path+orbit_file[:-3] +
                            f'_avg_x{fwhm_x}y{fwhm_y}z{fwhm_z}_nx{nx}.nc'),
                            mode="w", format="NETCDF4_CLASSIC")

#---------------------------_FIRST TEST ------------------------------------


#TEMP = TEMP[:,:,0:30].T
#x = data.x[0:30].values
#y = data.y.values
#z = data.height.values

# fwhm 10, 10, 10
#time_2 = time.time()
#avks_temp_dbl = apply_3d_kernel(TEMP, x, y, z, [10, 10, 10], only_kernel=False)
#elapsed_2 = time.time() - time_2
#np.save('avks_temp_dbl', avks_temp_dbl)

# fwhm 40, 10, 2
#time_3 = time.time()
#avks_temp_add = apply_3d_kernel(TEMP, x, y, z, [40, 10, 2], only_kernel=False)
#elapsed_3 = time.time() - time_3
#np.save('avks_temp_add', avks_temp_add)

# fwhm 5,5,5
#time_1 = time.time()
#avks_temp = apply_3d_kernel(TEMP, x, y, z, [2, 2, 2], only_kernel=False)
#elapsed_1 = time.time() - time_1
#np.save('avks_temp', avks_temp)

#print(time_1)
#print(time_2)
#print(time_3)

# %%

#fig, axs = plt.subplots(1,3)
#fig.suptitle('Vertically stacked subplots')
#axs[0].pcolormesh(TEMP[:,:,7])
#axs[1].pcolormesh(avks_temp_80_10_1_x15[7,:,:].T)
#img=axs[2].pcolormesh(TEMP[:,:,7]-avks_temp_80_10_1_x15[7,:,:].T)
#fig.colorbar(img, ax=axs[2])