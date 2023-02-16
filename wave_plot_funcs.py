# wave_plot_funcs
from datetime import datetime, timezone
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt

# hs time series plot
def plot_hs_ts(spotter_id, dt, hobs, hmod, show=True):
    rmse = np.sqrt( np.mean( (hobs - hmod)**2) )
    ts = 'RMSE: {:.2f} m'.format(rmse)
    fig, ax = plt.subplots(ncols=1)
    plt.plot(dt, hobs,linewidth=2,label='Measured')
    plt.plot(dt, hmod,linewidth=2,label='Modeled')
    plt.ylabel('Significant Wave Height $H_s$ (m)')
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
                                 horizontalalignment='right', verticalalignment='bottom')
    plt.legend(loc='upper left')
    plt.title(spotter_id)
    ax.tick_params(axis='x', rotation=45)
    plt.savefig(spotter_id+'_hs_ts.png', bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()
    

# hs scatter plot
def plot_hs_scat(spotter_id, hobs, hmod, show=True):
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
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
                 horizontalalignment='right', verticalalignment='bottom')
    plt.ylabel('Modeled $H_s$ (m)')
    plt.xlabel('Measured $H_s$ (m)')
    plt.axis([0, 12, 0, 12])
    plt.title(spotter_id)
    plt.savefig(spotter_id+'_hs_scat.png', bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
    
    
# tp time series plot
def plot_tp_ts(spotter_id, dt, tpobs, tpmod, show=True):
    rmse = np.sqrt( np.mean( (tpobs - tpmod)**2) )
    ts = 'RMSE: {:.2f} s'.format(rmse)
    fig, ax = plt.subplots(ncols=1)
    plt.plot(dt, tpobs, linewidth=2,label='Measured')
    plt.plot(dt, tpmod, linewidth=2,label='Modeled')
    plt.ylabel('Peak Wave Period $T_p$ (s)')
    plt.title(spotter_id)
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
                                 horizontalalignment='right', verticalalignment='bottom')
    plt.legend(loc='upper left')
    plt.ylim([2, 16])
    ax.tick_params(axis='x', rotation=45)
    plt.savefig(spotter_id+'_tp_ts.png', bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()


# tp scatter plot
def plot_tp_scat(spotter_id, tpobs, tpmod, show=True):
    gradient, intercept, r_value, p_value, std_err = stats.linregress(tpobs,tpmod)
    mn=np.min(tpobs)
    mx=np.max(tpobs)
    x1d=np.linspace(mn,mx,10)
    y1d=gradient*x1d+intercept
    ts = 'N={:.0f}\n$r^2$={:.2f}\nslope={:.2f}'.format( len(tpobs),r_value,gradient)
    fig, ax = plt.subplots(ncols=1)
    plt.plot([2,14],[2,14],'--',c='gray')
    ax.plot(x1d,y1d,'--',c='tab:blue')
    plt.plot(tpobs,tpmod,'.')
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
            horizontalalignment='right', verticalalignment='bottom')
    plt.ylabel('Modeled $T_p$ (s)')
    plt.xlabel('Measured $T_p$ (s)')
    plt.axis([2, 14, 2, 14])
    plt.title(spotter_id)
    plt.savefig(spotter_id+'_tp_scat.png', bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()

    
# dir time series plot
def plot_dir_ts(spotter_id, dt, mdirobs, mdirsobs, mdirmod, mdirsmod, show=True):
    rmse = np.sqrt( np.mean( (mdirobs - mdirmod)**2) )
    #axis = [-20,380,-20,380]
    ts = 'RMSE: {:.2f} $^\circ$'.format(rmse)
    fig, ax = plt.subplots(ncols=1)
    plt.errorbar(dt,mdirobs,yerr=mdirsobs,linewidth = 3,capsize=2, label='Measured')
    plt.errorbar(dt,mdirmod,yerr=mdirsmod,linewidth = 2, capsize=2, label='Modeled')
    plt.ylabel('Mean Wave Direction (${^\circ}$)')
    plt.title(spotter_id)
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
            horizontalalignment='right', verticalalignment='bottom')
    plt.legend(loc='upper left')
    plt.title(spotter_id)
    ax.tick_params(axis='x', rotation=45)
    plt.ylim([-20, 380])
    plt.savefig(spotter_id+'_dir_ts.png', bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()
    
# mdir (mean direction) scatter plot
def plot_mdir_scat(spotter_id, mdirobs, mdirmod, show=True):
    obs = mdirobs
    mod = mdirmod
    axis = [0,360,0,360]
    gradient, intercept, r_value, p_value, std_err = stats.linregress(obs,mod)
    mn=np.min(obs)
    mx=np.max(obs)
    x1d=np.linspace(mn,mx,10)
    y1d=gradient*x1d+intercept
    ts = 'N={:.0f}\n$r^2$={:.2f}\nslope={:.2f}'.format( len(obs),r_value,gradient)
    fig, ax = plt.subplots(ncols=1)
    plt.plot(axis[0:2],axis[0:2],'--',c='gray')
    ax.plot(x1d,y1d,'--',c='tab:blue')
    plt.plot(obs,mod,'.')
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
            horizontalalignment='right', verticalalignment='bottom')
    plt.ylabel('Modeled Mean Direction (${^\circ}$)')
    plt.xlabel('Measured Mean Direction (${^\circ}$)')
    plt.axis(axis)
    plt.title(spotter_id)
    plt.savefig(spotter_id+'_meandir_scat.png', bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()
    
# mdirs (spread) scatter plot
def plot_mdirs_scat(spotter_id, mdirsobs, mdirsmod, show=True):
    obs = mdirsobs
    mod = mdirsmod
    axis = [20,60,20,60]
    gradient, intercept, r_value, p_value, std_err = stats.linregress(obs,mod)
    mn=np.min(obs)
    mx=np.max(obs)
    x1d=np.linspace(mn,mx,10)
    y1d=gradient*x1d+intercept
    ts = 'N={:.0f}\n$r^2$={:.2f}\nslope={:.2f}'.format( len(obs),r_value,gradient)
    fig, ax = plt.subplots(ncols=1)
    plt.plot(axis[0:2],axis[0:2],'--',c='gray')
    ax.plot(x1d,y1d,'--',c='tab:blue')
    plt.plot(obs,mod,'.')
    plt.text(x=.95, y=.05, s=ts, fontsize=10, c='tab:blue',transform=ax.transAxes,\
            horizontalalignment='right', verticalalignment='bottom')
    plt.ylabel('Modeled Direc. Spread (${^\circ}$)')
    plt.xlabel('Measured Direc. Spread (${^\circ}$)')
    plt.axis(axis)
    plt.title(spotter_id)
    plt.savefig(spotter_id+'_mdirs_scat.png', bbox_inches='tight')
    if(show):
        plt.show()
    plt.close()