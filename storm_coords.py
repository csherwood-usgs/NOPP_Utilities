import numpy as np


def haversine(lon1, lat1, lon2, lat2, radius=6371.):
    """
    Calculate the great circle distance between two points
    on the earth. Input in decimal degrees.
    
    All args must be of equal length. 

    https://stackoverflow.com/questions/29545704/fast-haversine-approximation-python-pandas   
    Default radius of earth is  6371 km
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2  
    c = 2 * np.arcsin(np.sqrt(a))
    km = radius * c
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

def pcoord(x, y):
    """
    Convert x, y to polar coordinates r, az (geographic convention)
    r,az = pcoord(x, y)
    """
    r = np.sqrt(x**2 + y**2)
    az = np.degrees(np.arctan2(x, y))
    # az[where(az<0.)[0]] += 360.
    az = (az+360.)%360.
    return r, az

def xycoord(r, az):
    """
    Convert r, az [degrees, geographic convention] to rectangular coordinates
    x,y = xycoord(r, az)
    """
    x = r * np.sin(np.radians(az))
    y = r * np.cos(np.radians(az))
    return x, y

def dist_bearing(lon1, lat1, lon2, lat2):
    dist = haversine(lon1, lat1, lon2, lat2)
    brng = get_bearing(lon1, lat1, lon2, lat2)
    return dist, brng
