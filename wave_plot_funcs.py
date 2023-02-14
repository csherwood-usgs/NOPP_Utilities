# wave_plot_funcs
from datetime import datetime, timezone
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt

# hs time series plot
def plot_hs_ts(spotter_id, dt, hobs, hmod):
    rmse = np.sqrt( np.mean( (hobs - hmod)**2) )
    ts = 'RMSE: {:.2f} m'.format(rmse)
    fig, ax = plt.subplots(ncols=1)
    plt.plot(dt, hobs,linewidth=2,label='Measured')
    plt.plot(dt, hmod,linewidth=2,label='Modeled')
    plt.ylabel('Significant Wave Height $H_s$ (m)')
    plt.text(x=.05, y=.9, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes)
    plt.legend()
    ax.tick_params(axis='x', rotation=45)
    plt.savefig(spotter_id+'_hs_ts.png')
    plt.close()
    
# tp time series plot
def plot_tp_ts(spotter_id, dt, tpobs, tpmod):
    rmse = np.sqrt( np.mean( (tpobs - tpmod)**2) )
    ts = 'RMSE: {:.2f} s'.format(rmse)
    fig, ax = plt.subplots(ncols=1)
    plt.plot(dt, tpobs, linewidth=2,label='Measured')
    plt.plot(dt, tpmod, linewidth=2,label='Modeled')
    plt.ylabel('Peak Wave Period $T_p$ (s)')
    plt.title(spotter_id)
    plt.text(x=.05, y=.9, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes)
    plt.legend()
    ax.tick_params(axis='x', rotation=45)
    plt.savefig(spotter_id+'_tp_ts.png')
    #plt.show()
    plt.close()

# hs scatter plot
def plot_hs_scat(spotter_id, hobs, hmod):
    gradient, intercept, r_value, p_value, std_err = stats.linregress(hobs,hmod)
    mn=np.min(hobs)
    mx=np.max(hobs)
    x1d=np.linspace(mn,mx,10)
    y1d=gradient*x1d+intercept

    ts = 'N={:.0f}\n$r^2$={:.2f}\nslope={:.2f}'.format( len(hobs),r_value,gradient)

    fig, ax = plt.subplots(ncols=1)
    plt.plot([0,12],[0,12],'--',c='gray')
    ax.plot(x1d,y1d,'--',c='tab:blue')
    plt.plot(hobs,hmod,'.')
    plt.text(x=.05, y=.9, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes)
    plt.ylabel('Modeled $H_s$ (m)')
    plt.xlabel('Measured $H_s$ (m)')
    plt.axis([0, 12, 0, 12])
    plt.title(spotter_id)
    plt.savefig(spotter_id+'_hs_scat.png')
    plt.close()
    
# tp scatter plot
def plot_tp_scat(spotter_id, tpobs, tpmod):
    gradient, intercept, r_value, p_value, std_err = stats.linregress(tpobs,tpmod)
    mn=np.min(tpobs)
    mx=np.max(tpobs)
    x1d=np.linspace(mn,mx,10)
    y1d=gradient*x1d+intercept

    ts = 'N={:.0f}\n$r^2$={:.2f}\nslope={:.2f}'.format( len(tpobs),r_value,gradient)

    fig, ax = plt.subplots(ncols=1)
    plt.plot([4,14],[4,14],'--',c='gray')
    ax.plot(x1d,y1d,'--',c='tab:blue')
    plt.plot(tpobs,tpmod,'.')
    plt.text(x=.05, y=.7, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes)
    plt.ylabel('Modeled $T_p$ (s)')
    plt.xlabel('Measured $T_p$ (s)')
    plt.axis([4, 14, 4, 14])
    plt.title(spotter_id)
    plt.savefig(spotter_id+'_tp_scat.png')
    #plt.show()
    plt.close()