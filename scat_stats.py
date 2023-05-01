# statistics for scatter plots
import numpy as np

# Collection of statistics used to evaluate agreement between
# two sets of points (e.g. scatterplots).
# No attempt to optimize these has been performed.

# Might want to pass flattened arrays to these

# S = simulated (modeled) data points
# O = observed (measured, real) data points

# Statistics from MBCM
# Menaschi et al. (2013) Problems in RMSE-based wave model validations" Ocean Modeling 72:53-58.

def calc_HH( S, O ):
    """Equation MBCM 22
    The authors argue that this statistic is better than RMSE or NRMSE especially
    when null or negative biase and lots of scatter exist.
    Dimensionless.
    """
    return np.sqrt( np.nansum( ( S - O )**2 ) / np.nansum( S*O ) )

def calc_bias( S, O ):
    """Equation MBCM 6
    Units are same as input.
    """
    ok = np.isfinite(S+O) # eliminate if NaNs in either dataset
    return np.mean( S[ok] ) - np.mean( O[ok] )
           
def calc_RMSE( S, O ):
    """Equation MBCM  2
    Units are same as input.
    """
    return np.sqrt( np.nanmean( (S-O)**2 ) )

def calc_NRMSE( S, O ):
    """Equation MBCM 1
    Dimensionless.
    """
    ok = np.isfinite(S+O) # eliminate if NaNs in either dataset
    return np.sqrt( np.sum( (S[ok]-O[ok])**2 ) / np.sum( O[ok]**2 ) )

# Other stats
def calc_rho( S, O ):
    """Pearsons R
    Ranges from -1 (no correlation) to 1 (perfect correlation)
    Dimensionless.
    """
    ok = np.isfinite(S+O) # eliminate if NaNs in either dataset
    return  np.corrcoef( S[ok], O[ok] )[0,0]
    
def scat_stats_array( S, O ):
    S = S.flatten()
    O = O.flatten()
    assert S.shape == O.shape
    N = np.prod(O.shape)
    Nnan = np.sum(np.isnan(O))
    RMSE = calc_RMSE( S, O )
    rho = calc_rho( S, O )
    bias = calc_bias( S, O )
    NRMSE = calc_NRMSE( S, O )
    HH = calc_HH( S, O )
    return np.array([N, Nnan, RMSE, rho, bias, NRMSE, HH])

def scat_stats_string( S, O, sep_lines=True ):
    a = scat_stats_array( S, O )
    s = 'N: {0:.0f}\nNnan: {1:.0f}\nRMSE: {2:.3f}\n\rho: {3:.3f}\nBias: {4:.3f}\nNRMSE: {5:.3f}\nHH: {6:.3f}'.\
    format( a[0],a[1],a[2],a[3],a[4], a[5], a[6] ) 
    return s
    
    