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

def get_point_at_distance(lat1, lon1, d, bearing, R=6371):
    """
    lat: initial latitude, in degrees
    lon: initial longitude, in degrees
    d: target distance from initial (km)
    bearing: (true) heading in degrees
    R: optional radius of sphere, defaults to mean radius of earth

    Returns new lat/lon coordinate {d}km from initial, in degrees
    https://stackoverflow.com/questions/7222382/get-lat-long-given-current-point-distance-and-bearing
    """
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    a = np.radians(bearing)
    lat2 = np.arcsin(np.sin(lat1) * np.cos(d/R) + np.cos(lat1) * np.sin(d/R) * np.cos(a))
    lon2 = lon1 + np.arctan2(
        np.sin(a) * np.sin(d/R) * np.cos(lat1),
        np.cos(d/R) - np.sin(lat1) * np.sin(lat2)
    )
    return (np.degrees(lat2), np.degrees(lon2),)

def dist_bearing(lon1, lat1, lon2, lat2):
    dist = haversine(lon1, lat1, lon2, lat2)
    brng = get_bearing(lon1, lat1, lon2, lat2)
    return dist, brng
