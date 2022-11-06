from rawdata import time_tools as time_tools
import datetime as DT

def test_onboard_time():
    utctime = DT.datetime(2022,11,4,16,23,45,120000,tzinfo=DT.timezone.utc)
    onboard_time = time_tools.utc_to_onboardTime(utctime)
    
    assert onboard_time == 1351614243.12
    assert time_tools.onboardTime_to_utc(onboard_time) == utctime

if __name__ == "__main__":
    test_onboard_time()