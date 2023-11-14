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
    
def arc(ax, x0, y0, r, az0, az1, pstring='-', col='black', npts=20, alpha=.4, zorder=0,
        vmin=0., vmax=30., cmap=cm.Spectral_r ):
    """
    Draw an arc with radius r from az0 to az1 with npts
    """
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    az = np.linspace( az0, az1, npts )
    x, y = xycoord( r, az )
    xp, yp = x+x0, y+y0
    ax.plot(xp, yp, pstring, c=col, alpha=alpha, zorder=zorder )

    
def arcs(ax, x0, y0, r, az0, az1, pstring='-', col='black', npts=20, alpha=.4, zorder=0,
        vmin=0., vmax=30., cmap=cm.Spectral_r, rscale=True ):
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
    
def plt_spread(ax, x0, y0, f, S, dm, sprd, sf=2., ps = 25, fc=1./33., sfr=1., ec='black', vmin=-2.2, vmax=1.5, cmap='Reds', rscale=False, cbar=True):
    # Radial plot for a single observation, with or without colorbar.
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
    # cmap is colormap
    # rscale=True scales the spread according to the radius
    
    # for log scaling, suggest
    # vmin = -2.2
    # vmax = 1.5

    # for linear scaling, suggest
    # vmin=0.
    # vmax=30.
    
     # location of radial rings
    radii_f=[0.05, 0.1, 0.2, 0.3, 0.4]

    # sort the data so highest values will plot on top
    isort = np.argsort(S)
    S = S[isort]
    f = f[isort]
    dm = dm[isort]
    sprd = sprd[isort]

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(1, 1, layout='constrained')
    #ax = fig.add_subplot()
    ax.set_aspect('equal', adjustable='box')
    # frequency rings
    for fr in np.array( radii_f ):
        r = logr( fr )
        circle(ax, x0, y0, r, zorder=0)
    # label rings
    for i in np.array(radii_f):
        ptext(ax, x0, y0, logr(i), 45, "{}".format(i) )
    ax.axis('off')

    nf = np.shape(f)[0]
    sx = np.zeros_like(f)
    sy = np.zeros_like(f)
    for i in range(nf):
        sx[i], sy[i] = xycoord( logr( f[i] ), dm[i] )
        r = logr( f[i], fc=fc, sfr=sfr) 
        if rscale:
            sprd[i] = sprd[i]/r
        azlo = (dm[i]-sprd[i]/2.)
        azhi = (dm[i]+sprd[i]/2.)
        arc(ax, x0, y0, r , azlo, azhi, col=cmap(norm(S[i])),
            zorder=1, alpha = 0.7, vmin=vmin, vmax=vmax, cmap=cmap )

        ax.scatter( sx, sy, c = S, s = ps, vmin=vmin, vmax=vmax, edgecolors=ec, cmap=cmap, zorder=3, alpha = 0.8 )
    if(cbar):
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        ax.figure.colorbar(sm, ax = ax, shrink=.7, label = r'log$_{10}$[ Energy Density (m$^2$/Hz) ]')
        
def plt_2spread(ax, x0, y0, fobs, Sobs, dmobs, sprdobs, 
                            fmod, Smod, dmmod, sprdmod,
                sf=2., ps = 25, fc=1./33., sfr=1., ec='black', vmin=-2.2, vmax=1.5, cmap='Reds', rscale=False, cbar=True, titles=['Drifter','SWAN']):
    # Radial plot for a pair of model-data observations, with or without colorbar, with or without titles.
    # f = frequency [Hz]
    # S(f) = spectral density [m^2 (deg*Hz)^-1 ]
    # dm(f) = mean direction (deg)
    # sprd(f) = directional spread (deg)  
    # plot using radial coords
    # sf is scale factor for sigma
    # ps is point size for dot
    # fc is lowest frequency that can be plotted (center of circle)
    # sfr is general scaling factor
    # ec is edgecolor for dot
    # cmap is colormap
    # rscale=True scales the spread according to the radius
    
    # for log scaling, suggest
    # vmin = -2.2
    # vmax = 1.5

    # for linear scaling, suggest
    # vmin=0.
    # vmax=30.
    
     # location of radial rings
    radii_f=[0.05, 0.1, 0.2, 0.3, 0.4]

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(1, 2, layout='constrained')
    #ax = fig.add_subplot()
    for k, ax in enumerate(axs):
        ax.set_aspect('equal', adjustable='box')
        # frequency rings
        for fr in np.array( radii_f ):
            r = logr( fr )
            circle(ax, x0, y0, r, zorder=0)
        # label rings
        for i in np.array(radii_f):
            ptext(ax, x0, y0, logr(i), 45, "{}".format(i) )
        ax.axis('off')
        
        if k==0:
            # observations
            f = fobs
            S=Sobs
            dm=dmobs
            sprd=sprdobs
        else:
            # model
            f = fmod
            S=Smod
            dm=dmmod
            sprd=sprdmod
            
        # sort the data so highest values will plot on top
        isort = np.argsort(S)
        S = S[isort]
        f = f[isort]
        dm = dm[isort]
        sprd = sprd[isort]

        nf = np.shape(f)[0]
        sx = np.zeros_like(f)
        sy = np.zeros_like(f)
        for i in range(nf):
            sx[i], sy[i] = xycoord( logr( f[i] ), dm[i] )
            r = logr( f[i], fc=fc, sfr=sfr) 
            if rscale:
                sprd[i] = sprd[i]/r
            azlo = (dm[i]-sprd[i]/2.)
            azhi = (dm[i]+sprd[i]/2.)
            arc(ax, x0, y0, r , azlo, azhi, col=cmap(norm(S[i])),
                zorder=1, alpha = 0.7, vmin=vmin, vmax=vmax, cmap=cmap )

            ax.scatter( sx, sy, c = S, s = ps, vmin=vmin, vmax=vmax, edgecolors=ec, cmap=cmap, zorder=3, alpha = 0.8 )

            ax.set_title( titles[k] )
    if(cbar):
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        ax.figure.colorbar(sm, ax = axs[1], shrink=.4, label = r'log$_{10}$[ Energy Density (m$^2$/Hz) ]')
    
    
def plt_overlay_spread(ax, x0, y0, fobs, Sobs, dmobs, sprdobs, 
                            fmod, Smod, dmmod, sprdmod,
                sf=2., ps = 25, fc=1./33., sfr=1., ec='black', vmin=-2.2, vmax=1.5, cmap='Reds', rscale=False, cbar=True, title=['Drifter v. SWAN']):
    # Radial plot for a pair of model-data observations, with or without colorbar, with or without titles.
    # f = frequency [Hz]
    # S(f) = spectral density [m^2 (deg*Hz)^-1 ]
    # dm(f) = mean direction (deg)
    # sprd(f) = directional spread (deg)  
    # plot using radial coords
    # sf is scale factor for sigma
    # ps is point size for dot
    # fc is lowest frequency that can be plotted (center of circle)
    # sfr is general scaling factor
    # ec is edgecolor for dot
    # cmap is colormap
    # rscale=True scales the spread according to the radius
    
    # for log scaling, suggest
    # vmin = -2.2
    # vmax = 1.5

    # for linear scaling, suggest
    # vmin=0.
    # vmax=30.
    
     # location of radial rings
    radii_f=[0.05, 0.1, 0.2, 0.3, 0.4]

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(1, 1, layout='constrained')

    ax.set_aspect('equal', adjustable='box')
    # frequency rings
    for fr in np.array( radii_f ):
        r = logr( fr )
        circle(ax, x0, y0, r, zorder=0)
    # label rings
    for i in np.array(radii_f):
        ptext(ax, x0, y0, logr(i), 45, "{}".format(i) )
    ax.axis('off')
    for k in range(2):
        
        if k==0:
            # observations
            f = fobs
            S=Sobs
            dm=dmobs
            sprd=sprdobs
            cmap = cm.Reds
        else:
            # model
            f = fmod
            S=Smod
            dm=dmmod
            sprd=sprdmod
            cmap = cm.Greys
            
        # sort the data so highest values will plot on top
        isort = np.argsort(S)
        S = S[isort]
        f = f[isort]
        dm = dm[isort]
        sprd = sprd[isort]

        nf = np.shape(f)[0]
        sx = np.zeros_like(f)
        sy = np.zeros_like(f)
        for i in range(nf):
            sx[i], sy[i] = xycoord( logr( f[i] ), dm[i] )
            r = logr( f[i], fc=fc, sfr=sfr) 
            if rscale:
                sprd[i] = sprd[i]/r
            azlo = (dm[i]-sprd[i]/2.)
            azhi = (dm[i]+sprd[i]/2.)
            arc(ax, x0, y0, r , azlo, azhi, col=cmap(norm(S[i])),
                zorder=1, alpha = 0.7, vmin=vmin, vmax=vmax, cmap=cmap )
            ax.scatter( sx, sy, c = S, s = ps, vmin=vmin, vmax=vmax, edgecolors=ec, 
                       cmap=cmap, zorder=3, alpha = 0.8 )
            ax.set_title( title )
        if(cbar):
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            if k==0:
                ax.figure.colorbar(sm, ax = ax, shrink=.4, pad = .1,
                                   label = r'Drifter log$_{10}$[ Energy Density (m$^2$/Hz) ]')
            if k==1:
                ax.figure.colorbar(sm, ax = ax, shrink=.4, pad = .2,
                                   label = r'SWAN log$_{10}$[ Energy Density (m$^2$/Hz) ]')
    
    
def logr( f, fc =  1./33., sfr = 1. ):
    # calculations to convert frequency to plot radius
    # fc is lowest frequency that can be plotted (center of circle)
    # sfr is general scaling factor    
    r = sfr * (np.log10( f ) - np.log10( fc ))
    return r

def setup_radial_plot( cmap=cm.Spectral_r ):
    radii_f=[0.05, 0.1, 0.2, 0.3, 0.4]

    # placeholder in case we need to repostion the plot
    x0 = 0.
    y0 = 0.

    fig = plt.figure( )
    ax = fig.add_subplot()
    ax.set_aspect('equal', adjustable='box')

    # frequency rings
    for fr in np.array( radii_f ):
        r = logr( fr )
        circle(ax, x0, y0, r, zorder=0)

    # label rings
    for i in np.array(radii_f):
        ptext(ax, x0, y0, logr(i), 45, "{}".format(i) )
    ax.axis('off')

    cmap = cmap
    return fig, ax
