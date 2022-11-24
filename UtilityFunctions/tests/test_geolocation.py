from geolocation import satellite as satellite
import datetime as DT

def test_satellite():
    date = DT.datetime(2022,11,24,13,30,0)#        2022-11-23 08:20:00                  
    satlat,satlon,satLT,nadir_sza,nadir_mza,TPlat,TPlon,TPLT,TPsza,TPssa = satellite.get_position(date)

if __name__ == "__main__":
    test_satellite()