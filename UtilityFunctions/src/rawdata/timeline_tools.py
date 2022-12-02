import pandas as pd
import datetime as DT
import matplotlib.pyplot as plt


def load_schedule(filename):
    """Loads a timeline schedule overview and returns dataframe

    Arguments:
        filename (sting): name of .csv file to load

    Returns:
        df (:obj:`datetime`): dataframe with overview.

    """
        
    df = pd.read_csv(filename)
    df['start_date'] = pd.to_datetime(df['start_date'],tzinfo=DT.timezone.utc)
    df['end_date'] = pd.to_datetime(df['end_date'],tzinfo=DT.timezone.utc)

    return df

def plot_schedule(df,column='name'):
    """Plots info from timeline schedule

    Arguments:
        df (obj:`dataframe`): Pandas dataframe holding the schedule
        column: Value to plot (default: name)
    Returns:
        None

    """
    for x1, x2, y in zip(df["start_date"], df["end_date"], df[column]):
        plt.plot([x1, x2], [y, y])
    