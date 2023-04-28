# statistics for scatter plots
import numpy as np

# might want to pass flattened arrays to these

# S = simulated (modeled) data points
# O = observed (measured, real) data points

def calc_HH( S, O ):
    """Equation 22"""
    return np.sqrt( np.nansum( ( S - O )**2 ) / np.nansum( S*O ) )

def calc_bias( S, O ):
    """Equation 6"""
    ok = np.isfinite(S+O) # eliminate if NaNs in either dataset
    return np.mean( S[ok] ) - np.mean( O[ok] )
    