"""
Codes for converting among various directional wave formats and calculating bulk statistics.

Many of these routines include code and/or guidance from Isabel Houghton at Sofar Ocean.

csherwood@usgs.gov
"""
import numpy as np

def integrate_in_frequency(data_field, frequencies, faxis=0):
    return np.trapz(data_field, frequencies, axis=faxis).squeeze()


def integrate_in_direction(data_field, directional_bin_width, daxis=1):
    return np.sum(data_field * directional_bin_width, axis=daxis).squeeze()


def integrate_frequency_and_direction(data_field, frequencies, directional_bin_width, faxis=0, daxis=1):
    fsum = integrate_in_frequency(data_field, frequencies, faxis=faxis).squeeze()
    return integrate_in_direction(fsum, directional_bin_width, daxis=daxis).squeeze()


def to_moment( R1,R2, alpha1,alpha2 ):
    """
    Convert NDBC directional parameters r1, r2, theta1, theta2 to Fourier coeffs. a1, a2, b1, and b2.
    
    This code was written by Pieter Smit and provided courtesy of Isabel Houghton, Sofar Ocean.
    
    Note: from the NDBC website
   
    The R1 and R2 values in the monthly and yearly historical data files are scaled by 100,
    a carryover from how the data are transported to the archive centers. The units are
    hundredths, so the R1 and R2 values in those files should be multiplied by 0.01.
   
    May have to multiply R1/R2 with a 100 depending on source
   
    :param R1:
    :param R2:
    :param alpha1:
    :param alpha2:
    :return:
    """
    to_rad = np.pi / 180
    eps = 1e-16
    angle1 = (270 - alpha1) * to_rad
    angle2 = (540 - 2 * alpha2) * to_rad

    a0 = 1/(2*np.pi) + eps
    a1 =R1 * np.cos(angle1) / a0
    b1 =R1 * np.sin(angle1) / a0
    a2 =R2 * np.cos(angle2) / a0
    b2 =R2 * np.sin(angle2) / a0

    return a1, b1, a2, b2


def to_Fourier( Efth, frequencies, direction_radians, directional_bin_width_deg, faxis=0, daxis=1 ):
    """
    Convert 2d spectrum to Fourier coefficients
    Based on code from Isabel Houghton, Sofar
    """
    eps = 1e-16
    Ef = integrate_in_direction(Efth, directional_bin_width_deg, daxis=daxis) + eps
    print('shape of Ef: ',np.shape(Efth))
    print('shape of Ef: ',np.shape(Ef))

    # each directional coefficient (a1, b1..) is a Fourier coefficient of the 2D WW3 spectrum,
    # normalized by the 1D (directionally-integrated) energy density spectrum
    a1 = integrate_in_direction(Efth * np.cos(direction_radians),
                                directional_bin_width=directional_bin_width_deg, daxis=daxis) / Ef

    a2 = integrate_in_direction(Efth * np.cos(2 * direction_radians),
                                directional_bin_width=directional_bin_width_deg, daxis=daxis) / Ef

    b1 = integrate_in_direction(Efth * np.sin(direction_radians),
                                directional_bin_width=directional_bin_width_deg, daxis=daxis) / Ef

    b2 = integrate_in_direction(Efth * np.sin(2 * direction_radians),
                                directional_bin_width=directional_bin_width_deg, daxis=daxis) / Ef
    
    return Ef, a1, a2, b1, b2


#### 1d statistics from 1d data

def calc_m0_1d( spec1d, frequencies ):
    return np.trapz( spec1d, frequencies )


def calc_m1_1d( spec1d, frequencies ):
    return np.trapz( spec1d * frequencies, frequencies )


def calc_Hs_1d( spec1d, frequencies ):
    return 4. * np.sqrt( calc_m0_1d( spec1d, frequencies ) )


def calc_spec1d_2d( spec2d, directional_bin_width, daxis=1 ):
    # integrate over direction to get 1d spec from 2d spec
    return integrate_in_direction( spec2d, directional_bin_width, daxis=daxis ) 


def calc_Hs_2d( data_field, frequencies, directional_bin_width, faxis=0, daxis=1 ):
    return( 4.*np.sqrt( integrate_frequency_and_direction(data_field, frequencies, directional_bin_width, faxis=faxis, daxis=daxis).squeeze()) )


def calc_sigmaf_1d( spec1d, frequencies ):
    # frequency spread () LeMerle et al., 2021 eqn. 3 based on Blackman & Tukey (1959)
    return np.trapz( spec1d, frequencies )**2 / np.trap1d( spec1d**2, frequencies )


def calc_Qp_1d( spec1d, frequencies ):
    # frequency peakedness () Le Merle et al., 2021 eqn 4 based on Goda (1976)
    return 2.*np.trapz( spec1d**2, frequencies ) / (np.trapz( spec1d, frequencies ))**2 


def calc_theta_mean_a1b1( a1, b1, frequencies ):
    # mean direction (degrees)
    # no weighting
    # Le Merle et al., 2021 eqn 5; Isabels code, SoFar manual
    a1bar = np.trapz(a1, frequencies)
    b1bar = np.trapz(b1, frequencies)
    thetam = 270. - (180./np.pi) * np.arctan2( b1bar, a1bar )
    return thetam


def calc_weighted_theta_mean_a1b1( a1, b1, spec1d, frequencies ):
    # mean direction (degrees)
    # no weighting
    # Le Merle et al., 2021 eqn 5; Isabels code, SoFar manual
    a1bar = np.trapz(a1 * spec1d, frequencies) / np.trapz( spec1d, frequencies )
    b1bar = np.trapz(b1 * spec1d,  frequencies) / np.trapz( spec1d, frequencies )
    thetamw = 270. - (180./np.pi) * np.arctan2( b1bar, a1bar )
    return thetamw


def calc_sigma_theta_a1b1( a1, b1, spec1d, frequencies ):
    # directional spread (degrees)
    # Herbers et al., eqn 3a; Le Merle et al., 2021 eqns. 5 & 6
    a1bar = np.trapz(a1 * spec1d, frequencies) / np.trapz( spec1d, frequencies )
    b1bar = np.trapz(b1 * spec1d,  frequencies) / np.trapz( spec1d, frequencies )
    s1 = (180./np.pi)*np.sqrt( 2. * ( 1.-np.sqrt( a1bar**2 + b1bar**2 ) ) )
    s1[s1>=360.]=s1[s1>=360.]-360.
    s1[s1<0.]=s1[s1<0.]+360.
    return s1


def calc_sigma_theta2_a2b2( a2, b2, spec1d, frequencies ):
    # second estimate of directional spread (f) (degrees)
    # Herbers et al., eqn 3a; Ardhuin et al., eqns 5 & 6.
    a2bar = np.trapz(a2 * spec1d, frequencies) / np.trapz( spec1d, frequencies )
    b2bar = np.trapz(b2 * spec1d,  frequencies) / np.trapz( spec1d, frequencies )
    s2 = (180./np.pi)*np.sqrt( 0.5 * ( 1.-np.sqrt( a2bar**2 + b2bar**2 ) ) )
    s2[s2>=360.]=s2[s2>=360.]-360.
    s2[s2<0.]=s2[s2<0.]+360.
    return s2


def calc_FSPR_1d( spec1d, frequencies, TM02 ):
    # SWAN Manual; Battjes and Van Vledder (1984), their eqn. 5
    i = 0. + 1.j
    omega = 2.*np.pi*frequencies
    return np.abs( np.trapz( spec1d*np.exp(i*omega*TM02), frequencies ).squeeze() ) / np.trapz(spec1d, frequencies ).squeeze()


def calc_TM01_2d( data_field, frequencies, directional_bin_width, faxis=0, daxis=1 ):
    dsum = integrate_in_direction(data_field, directional_bin_width, daxis=daxis).squeeze()
    return ( integrate_in_frequency(dsum, frequencies, faxis=0).squeeze()/integrate_in_frequency(frequencies*dsum, frequencies, faxis=0).squeeze() )


def calc_TM02_2d( data_field, frequencies, directional_bin_width, faxis=0, daxis=1 ):
    dsum = integrate_in_direction(data_field, directional_bin_width, daxis=daxis).squeeze()
    return np.sqrt( integrate_in_frequency(dsum, frequencies, faxis=0).squeeze()/integrate_in_frequency(frequencies*frequencies*dsum, frequencies, faxis=0).squeeze() )


def calc_TM01_1d( spec1d, frequencies, faxis=0 ):
    return np.traz( spec1d, frequencies, axis=faxis ).squeeze() / np.trapz( frequencies*spec1d, frequencies, axis=faxis).squeeze()


def calc_TM02_1d( spec1d, frequencies, faxis=0 ):
    return np.sqrt( np.traz( spec1d, frequencies, axis=faxis ).squeeze() / np.trapz( frequencies*frequencies*spec1d, frequencies, axis=faxis).squeeze() )


#def calc_DSPR_spec2d( spec2d, 