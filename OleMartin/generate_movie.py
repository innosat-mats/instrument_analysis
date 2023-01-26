#%% Import modules
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot

#%% Select on explicit time
start_time = DT.datetime(2023,1,8,19,0,0)
stop_time = DT.datetime(2023,1,8,19,5,0)

#%% Get MATS data as a pandas dataframe
#df = read_MATS_data(start_time,stop_time)


#%% plot a nice image
#plot_image(df.iloc[3])

#%% plot_all_images and write to file using Bjorns plotting function
#simple_plot(df,'./',custom_cbar=True,ranges=[500,20000])

#%% plot_all_images with orbit data and write to file using Bjorns plotting function
import subprocess

def run_ffmpeg(infiles,outfile,framerate=30):   
    ffmpeg = "ffmpeg"

    commands = ' '.join([ffmpeg,
        "-r ", str(framerate),
        "-pattern_type glob ",
        "-i ", infiles,
        "-c:v ", "libx264 ",
        "-preset ",
        "-pix_fmt ",
        "yuv420p ",
        outfile
        ])
   
    print(commands)
    if subprocess.run(commands).returncode == 0:
        print ("FFmpeg Script Ran Successfully")
    else:
        print ("There was an error running your FFmpeg script")


outdir = './'
#orbit_plot(df,outdir)
run_ffmpeg(outdir + 'CCDSEL1/*.png','out')