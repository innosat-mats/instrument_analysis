#%%
from rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from plotting.plotCCD import simple_plot

schedule = pd.read_csv('/home/olemar/Projects/Universitetet/MATS/instrument_analysis/Operational_analysis/commisioning/timeline_schedule.csv')

start_time = DT.datetime.strptime(schedule[schedule.name=="ODIN1"].start_date.values[0],'%Y-%m-%d %H:%M:%S')
stop_time = DT.datetime.strptime(schedule[schedule.name=="ODIN1"].end_date.values[0],'%Y-%m-%d %H:%M:%S')
#%%
df = read_MATS_data(start_time,stop_time)
simple_plot(df,'./',custom_cbar=True,ranges=[300,1000]) #plot_all_images
#%%


