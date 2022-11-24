import datetime as DT


def select_on_time(start_date,end_date,df):
    """Function filters a MATS dataframe based on start and end time

    Arguments:
        start_date (:obj:`datetime`):
        end_date (:obj:`datetime`):
        df (:obj:`pandas dataframe`):

    Returns:
        df (:obj:`pandas dataframe`): filtered MATS dataframe.

    """

    if start_date.tzinfo == None:
        start_date = start_date.replace(tzinfo=DT.timezone.utc)
    if end_date.tzinfo == None:
        end_date = end_date.replace(tzinfo=DT.timezone.utc)

    df = df[df.TMHeaderTime>start_date]
    df = df[df.TMHeaderTime<end_date]

    return df