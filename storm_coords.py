import numpy as np


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    
    All args must be of equal length. 

    https://stackoverflow.com/questions/29545704/fast-haversine-approximation-python-pandas   
    Changed radius of earth to 6371. CRS
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6371. * c
    return km

def get_bearing(lon1, lat1, lon2, lat2):
    #https://stackoverflow.com/questions/54873868/python-calculate-bearing-between-two-lat-long
    dLon = (lon2 - lon1)
    x = np.cos(np.radians(lat2)) * np.sin(np.radians(dLon))
    y = np.cos(np.radians(lat1)) * np.sin(np.radians(lat2)) - np.sin(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(dLon))
    brng_rad = np.arctan2(x,y)
    brng_deg = np.degrees(brng_rad)

    return brng_deg

def dist_bearing(lon1, lat1, lon2, lat2):
    dist = haversine(lon1, lat1, lon2, lat2)
    brng = get_bearing(lon1, lat1, lon2, lat2)
    return dist, brng

if __name__ == '__main__':
    lat1 = 47.1
    lon1 = -62.
    lat2 = 47.2
    lon2 = 62.1

    d, az = dist_bearing( lon1, lat1, lon2, lat2)
    print(d, az)