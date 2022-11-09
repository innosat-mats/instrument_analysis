#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:53:36 2022

@author: lindamegner


This is more or less a duplicate of read_and_calibrate_all_files_in_directory but it is meant to be used as a script.
"""

import matplotlib.pyplot as plt
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.get_temperature import create_temperature_info_array
#import rawdata.time_tools as time_tools
from mats_l1_processing.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems

def plot_full_temperature_info(temperaturedata, relativetimedata):
    HTR1A = temperaturedata[:, 0]
    HTR1B = temperaturedata[:, 1]
    HTR2A = temperaturedata[:, 2]
    HTR2B = temperaturedata[:, 3]
    HTR8A = temperaturedata[:, 4]
    HTR8B = temperaturedata[:, 5]

    plt.plot(relativetimedata / 60.0, HTR1A,
             label="splitter plate, regulation")
    plt.plot(relativetimedata / 60.0, HTR1B, label="splitter plate, measuring")
    plt.plot(relativetimedata / 60.0, HTR2A, label="limb house, regulation")
    plt.plot(relativetimedata / 60.0, HTR2B, label="limb house, measuring")
    plt.plot(relativetimedata / 60.0, HTR8A, label="UV2 CCDn")
    plt.plot(relativetimedata / 60.0, HTR8B, label="UV1 CCDn")
    plt.xlabel("Time since epoch [min]")
    plt.ylabel("Temperature [C]")
    plt.legend()
    plt.show()
    plt.savefig("HTRmeasurements.jpg")




#directory='/Users/lindamegner/MATS/retrieval/FlightData/221106_payload_first_switchon/RacFiles_out'
directory='/Users/lindamegner/MATS/retrieval/FlightData/221109_first_images/RacFiles_out/'
#directory = '/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/testdata/binning_test_20200812/RacFiles_out/'
#add_rac_temp_data(directory + "/HTR.csv", CCDitem, labtemp=999)
temperaturedata, relativetimedata = create_temperature_info_array(directory + "/HTR.csv")

# utctime = []
# for i in range(len(relativetimedata)):
#     utctime.append(time_tools.onboardTime_to_utc(relativetimedata[i]))

#plot_full_temperature_info(temperaturedata,utctime)

plot_full_temperature_info(temperaturedata, relativetimedata)

CCDitems = read_CCDitems(directory)
if calibrate:
    calibrate_CCDitems(CCDitems, instrument)

fig, ax = plt.subplots(1)
plotCCDitem(CCDitems[0], fig, ax, title='Fist Image !!!')