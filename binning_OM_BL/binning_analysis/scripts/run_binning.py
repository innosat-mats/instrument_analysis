"""
Created on Fri Aug 14 15:15:36 2020

@author: olemar
"""
# %%
# fmt: off
import sys
sys.path.insert(0, "/home/olemar/Projects/MATS/instrument_analysis/binning_OM_BL/binning_analysis/src/binning_tools/")

import binning_functions as bf
import pandas as pd
from matplotlib import pyplot as plt
from mats_l1_processing import read_in_functions
import numpy as np
from scipy.stats import binned_statistic
from matplotlib.pyplot import cm

# fmt: on


def filter_on_time(CCDitems, starttime=None, stoptime=None):
    I = []
    for i in range(len(CCDitems)):
        image_time = pd.to_datetime(
            CCDitems[i]["EXP Date"], format="%Y-%m-%dT%H:%M:%S.%fZ"
        )
        if (starttime != None) and (stoptime != None):
            if (image_time > starttime) and (image_time < endtime):
                I.append(i)
        elif (starttime != None) and (stoptime == None):
            if image_time > starttime:
                I.append(i)
        elif (starttime == None) and (stoptime != None):
            if image_time < starttime:
                I.append(i)
        else:
            Warning("Start or end time invalid")

    CCDitems = [CCDitems[i] for i in I]
    return CCDitems


def fit_curve(man_tot, inst_tot, deg, threshold=np.inf):

    all_simulated = man_tot.flatten()
    all_measured = inst_tot.flatten()

    # fit linear part
    low_simulated = all_simulated[all_simulated < threshold]
    low_measured = all_measured[all_simulated < threshold]

    low_measured_mean, bin_edges = binned_statistic(
        low_simulated, low_measured, "median", bins=2000
    )[0:2]

    bin_center = (bin_edges[1:] + bin_edges[0:-1]) / 2
    p_low = np.polyfit(
        bin_center[~np.isnan(low_measured_mean)],
        low_measured_mean[~np.isnan(low_measured_mean)],
        deg,
        full=True,
    )[0]

    return p_low, bin_center, low_measured_mean


def get_linearity(CCDitems, testtype, plot=True):
    testtype = "col"
    channels = [1, 2, 3, 4, 5, 6]
    plotting_factor = 10
    threshold = 4e3
    deg = 1

    color = cm.rainbow(np.linspace(0, 1, 7))

    for i in range(len(channels)):
        (
            man_tot,
            inst_tot,
            channel,
            test_type,
        ) = bf.get_binning_test_data_from_CCD_item(
            CCDitems, test_type_filter=testtype, channels=[channels[i]]
        )

        p_low, bin_center, low_measured_mean = fit_curve(
            man_tot, inst_tot, deg, threshold
        )

        if plot:
            plt.plot(
                man_tot.flatten()[::plotting_factor],
                inst_tot[::plotting_factor],
                ".",
                alpha=0.1,
                markeredgecolor="none",
                c=color[i],
            )
            plt.plot(
                np.arange(0, 40000),
                np.polyval(p_low, np.arange(0, 40000)),
                "-",
                c=color[i],
            )

        print(p_low)
    if plot:
        plt.show()

    return p_low


# %%

dirname = "/home/olemar/Projects/MATS/MATS-data/binning_test_20200812_racfiles/binning/"

CCDitems = read_in_functions.read_CCDitems(dirname)
print(len(CCDitems))

starttime = pd.to_datetime("2020-08-10T12:42Z", format="%Y-%m-%dT%H:%MZ")
endtime = pd.to_datetime("2020-08-16T13:42Z", format="%Y-%m-%dT%H:%MZ")

# CCDitems = filter_on_time(CCDitems, starttime, endtime)

get_linearity(CCDitems, "col", plot=True)

# %%
