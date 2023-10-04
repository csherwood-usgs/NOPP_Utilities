import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# spec_plot_funcs - Routines to facilitate plotting wave spectra

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

def circle(ax, x0, y0, r, pstring='--', col='gray', alpha=.6, npts=100, zorder=0):
    """
    Draw a circle centered at x0, y0 with radius r and npts
    """
    az = np.linspace(0., 360.)
    x, y = xycoord(r, az)
    xp, yp = x+x0, y+y0
    ax.plot(xp, yp, pstring, c=col, alpha=alpha, zorder=zorder)
    
def arc(ax, x0, y0, r, az0, az1, S, pstring='-', col='black', npts=20, alpha=.4, zorder=0,
        vmin=0., vmax=30., cmap=cm.Spectral_r ):
    """
    Draw an arc with radius r from az0 to az1 with npts
    """
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    az = np.linspace( az0, az1, npts )
    x, y = xycoord( r, az )
    xp, yp = x+x0, y+y0
    ax.plot(xp, yp, pstring, c=col, alpha=alpha, zorder=zorder )
    
def pline(ax, x0, y0, x, y, pstring='-', col='black', alpha=.6, zorder=1):
    """
    Draw a line using arrays x, y
    """
    ax.plot( x0+x, x0+y, pstring, c=col, alpha=alpha, zorder=zorder)
    
def ptext(ax, x0, y0, r, az, txt, alpha=.6, zorder=0, fontsize=8 ):
    """
    Draw text using polar coordinates
    """
    x, y = xycoord( r, az )
    ax.text(x+x0, y+y0, txt, alpha=alpha, fontsize=fontsize )
    
def plt_data(ax, lon, lat, Hs, Tp, sigf, Dp, Dsprd):
    sf = 10 # scale factor for sigmaf
    ps = 50 # point size for Hs
    # coordinates for plot, based on peak period and peak direction
    wx, wy = xycoord(Tp, Dp)

    # coordinates for frequency spread, scaled by sf
    rlo, rhi = Tp-sigf*sf, Tp+sigf*sf
    rlox, rloy = xycoord(rlo, Dp)
    rhix, rhiy = xycoord(rhi, Dp)
    rx = np.array((rlox, rhix))
    ry = np.array((rloy, rhiy))
    arc(ax, x0, y0, Tp, Dp-Dsprd/2., Dp+Dsprd/2., zorder=1)
    pline(ax, x0, y0, rx, ry, zorder=1 )
    ax.scatter(wx, wy, ps, Hs, vmin=0., vmax=10., zorder=3)
    
def plt_rdata(ax, lon, lat, Hs, Tp, sigf, Dp, Dsprd, sf=2., ps = 85, fc=1./33., sfr=1., ec='black'):
    # plot using radial coords
    # sf is scale factor for sigma
    # ps is point size for dot
    # fc is lowest frequency that can be plotted (center of circle)
    # sfr is general scaling factor
    # ec is edgecolor for dot

    # coordinates for plot, based on peak period and peak direction
    wx, wy = xycoord(logr( 1./Tp, fc=fc, sfr=sfr), Dp)

    # coordinates for frequency spread, scaled by sf
    delT = sigf*Tp*sf
    # print(sigf, delT, Tp-delT, Tp+delT)
    flo = 1./(Tp - delT)
    fhi = 1./(Tp + delT)
    # print(flo, fhi)
    rlo, rhi = logr( flo, fc=fc, sfr=sfr ), logr( fhi, fc=fc, sfr=sfr )
    rlox, rloy = xycoord(rlo, Dp)
    rhix, rhiy = xycoord(rhi, Dp)
    rx = np.array((rlox, rhix))
    ry = np.array((rloy, rhiy))
    arc(ax, x0, y0, logr( 1./Tp, fc=fc, sfr=sfr) , Dp-Dsprd/2., Dp+Dsprd/2., zorder=1)
    pline(ax, x0, y0, rx, ry, zorder=1 )
    ax.scatter(wx, wy, ps, Hs, vmin=0., vmax=10., zorder=3, edgecolor=ec)
    
def plt_spread(ax, lon, lat, f, S, dm, sprd, sf=2., ps = 25, fc=1./33., sfr=1., ec='black', cmap='Reds'):
    # f = frequency [Hz]
    # S(f) = spectral density [m^2 (deg*Hz)^-1 ]
    # dirm(f) = mean direction (deg)
    # sprd(f) = directional spread (deg)
    
    # plot using radial coords
    # sf is scale factor for sigma
    # ps is point size for dot
    # fc is lowest frequency that can be plotted (center of circle)
    # sfr is general scaling factor
    # ec is edgecolor for dot
    vmin=0.
    vmax=30.
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    
    nf = np.shape(f)[0]
    sx = np.zeros_like(f)
    sy = np.zeros_like(f)
    for i in range(nf):
        sx[i], sy[i] = xycoord( logr( f[i] ), dm[i] )
        arc(ax, x0, y0, logr( f[i], fc=fc, sfr=sfr) , dm[i]-sprd[i]/2., dm[i]+sprd[i]/2., S,
            zorder=1, alpha = 0.6, vmin=vmin, vmax=vmax, cmap=cmap )
        
    ax.scatter( sx, sy, c = S, s = ps, vmin=vmin, vmax=vmax, cmap=cmap, zorder=3, alpha = 0.8 )

    
def logr( f, fc =  1./33., sfr = 1. ):
    # calculations to convert frequency to plot radius
    # fc is lowest frequency that can be plotted (center of circle)
    # sfr is general scaling factor    
    r = sfr * (np.log10( f ) - np.log10( fc ))
    return r
