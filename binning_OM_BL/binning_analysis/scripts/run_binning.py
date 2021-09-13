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

# fmt: on


# %%


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


dirname = "/home/olemar/Projects/MATS/MATS-data/BinningFlatfieldsIR3_100920/BinningIR3/"

CCDitems = read_in_functions.read_CCDitems(dirname)
print(len(CCDitems))

starttime = pd.to_datetime("2021-09-10T12:42Z", format="%Y-%m-%dT%H:%MZ")
endtime = pd.to_datetime("2021-09-10T13:42Z", format="%Y-%m-%dT%H:%MZ")

CCDitems = filter_on_time(CCDitems, starttime, endtime)
# %%

print(len(CCDitems))

(
    man_tot_exp,
    inst_tot_exp,
    channel_tot,
    test_type_tot,
) = bf.get_binning_test_data_from_CCD_item(CCDitems, test_type_filter="col")

# man_tot_col, inst_tot_col, channel_tot, test_type_tot = bf.get_binning_test_data(
#     dirname, test_type_filter="row", channels=[3]
# )

# man_tot_row, inst_tot_row, channel_tot, test_type_tot = bf.get_binning_test_data(
#     dirname, test_type_filter="col", channels=[3]
# )

# plt.plot(
#     man_tot_row.flatten(),
#     inst_tot_row / man_tot_row.flatten(),
#     ".",
#     alpha=0.01,
#     markeredgecolor="none",
# )

plt.plot(
    man_tot_exp.flatten(),
    inst_tot_exp,
    ".",
    alpha=0.01,
    markeredgecolor="none",
)
# plt.plot(
#     man_tot_col.flatten(),
#     inst_tot_col / man_tot_col.flatten(),
#     ".",
#     alpha=0.01,
#     markeredgecolor="none",
# )

plt.show()

print(test_type_tot)


# %%
