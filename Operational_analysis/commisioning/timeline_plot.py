# %%
import rawdata.timeline_tools as timeline_tools
import datetime as DT

# %%
directory='/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/commisioning/'
df = timeline_tools.load_schedule(directory + 'timeline_schedule.csv')
timeline_tools.plot_schedule(df)
# %%
