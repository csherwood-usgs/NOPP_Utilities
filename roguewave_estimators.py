# roguewave_estimators
# This code was copied from the Sofar Roguewave github site because I could not install that package on Win64.

### estimators.py
import numpy
#from roguewave.wavespectra.estimators.mem2 import mem2
#from roguewave.wavespectra.estimators.mem import mem
from numba_progress import ProgressBar
from typing import Literal

Estimators = Literal["mem", "mem2"]


# -----------------------------------------------------------------------------
#                       Boilerplate Interfaces
# -----------------------------------------------------------------------------
def estimate_directional_spectrum_from_moments(
    e: numpy.ndarray,
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
    direction: numpy.ndarray,
    method: Estimators = "mem2",
    **kwargs,
) -> numpy.ndarray:
    """
    Construct a 2D directional distribution based on the directional moments and a spectral
    reconstruction method.

    :param number_of_directions: length of the directional vector for the
    2D spectrum. Directions returned are in degrees

    :param method: Choose a method in ['mem','mem2']
        mem: maximum entrophy (in the Boltzmann sense) method
        Lygre, A., & Krogstad, H. E. (1986). Explicit expression and
        fast but tends to create narrow spectra anderroneous secondary peaks.

        mem2: use entrophy (in the Shannon sense) to maximize. Likely
        best method see- Benoit, M. (1993).

    REFERENCES:
    Benoit, M. (1993). Practical comparative performance survey of methods
        used for estimating directional wave spectra from heave-pitch-roll data.
        In Coastal Engineering 1992 (pp. 62-75).

    Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the
        directional distribution in ocean wave spectra.
        Journal of Physical Oceanography, 16(12), 2052-2060.

    """
    return (
        estimate_directional_distribution(a1, b1, a2, b2, direction, method, **kwargs)
        * e[..., None]
    )


def estimate_directional_distribution(
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
    direction: numpy.ndarray,
    method: Estimators = "mem2",
    **kwargs,
) -> numpy.ndarray:
    """
    Construct a 2D directional distribution based on the directional moments and a spectral
    reconstruction method.

    :param number_of_directions: length of the directional vector for the
    2D spectrum. Directions returned are in degrees

    :param method: Choose a method in ['mem','mem2']
        mem: maximum entrophy (in the Boltzmann sense) method
        Lygre, A., & Krogstad, H. E. (1986). Explicit expression and
        fast but tends to create narrow spectra anderroneous secondary peaks.

        mem2: use entrophy (in the Shannon sense) to maximize. Likely
        best method see- Benoit, M. (1993).

    REFERENCES:
    Benoit, M. (1993). Practical comparative performance survey of methods
        used for estimating directional wave spectra from heave-pitch-roll data.
        In Coastal Engineering 1992 (pp. 62-75).

    Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the
        directional distribution in ocean wave spectra.
        Journal of Physical Oceanography, 16(12), 2052-2060.

    """

    # Jacobian to transform distribution as function of radian angles into
    # degrees.
    Jacobian = numpy.pi / 180
    direction_radians = direction * Jacobian

    if method.lower() in ["maximum_entropy_method", "mem"]:
        # reconstruct the directional distribution using the maximum entropy
        # method.
        function = mem
    elif method.lower() in ["maximum_entrophy_method2", "mem2"]:
        function = mem2
    else:
        raise Exception(f"unsupported spectral estimator method: {method}")

    output_shape = list(a1.shape) + [len(direction)]
    if a1.ndim == 1:
        input_shape = [1, a1.shape[-1]]
    else:
        input_shape = [numpy.prod(a1.shape[0:-1]), a1.shape[-1]]

    a1 = a1.reshape(input_shape)
    b1 = b1.reshape(input_shape)
    a2 = a2.reshape(input_shape)
    b2 = b2.reshape(input_shape)

    number_of_iterations = a1.shape[0]
    if number_of_iterations < 10:
        disable = True
    else:
        disable = False

    if method != "mem2":
        msg = f"Reconstructing 2d spectrum with {method} using implementation: "
    else:
        solution_method = kwargs.get("solution_method", "scipy")
        msg = f"Reconstructing 2d spectrum with {method} using solution_method {solution_method}"

    with ProgressBar(total=number_of_iterations, disable=disable, desc=msg) as progress:
        res = function(direction_radians, a1, b1, a2, b2, progress, **kwargs)

    return res.reshape(output_shape) * Jacobian


### logliklihood.py
import numpy
from qpsolvers import solve_ls
#from .utils import get_direction_increment, get_rhs, get_constraint_matrix


def log_likelyhood(
    directions_radians: numpy.ndarray,
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
    progress,
    **kwargs
) -> numpy.ndarray:
    number_of_frequencies = a1.shape[-1]
    number_of_points = a1.shape[0]

    directional_distribution = numpy.zeros(
        (number_of_points, number_of_frequencies, len(directions_radians))
    )

    for ipoint in range(0, number_of_points):
        progress.update(1)
        directional_distribution[ipoint, :, :] = _log_likelyhood(
            directions_radians,
            a1[
                ipoint,
                :,
            ],
            b1[ipoint, :],
            a2[ipoint, :],
            b2[ipoint, :],
        )
    return directional_distribution


def _log_likelyhood(directions_radians: numpy.ndarray, a1, b1, a2, b2) -> numpy.ndarray:
    """
    Return the directional distribution that minimizes the variance (D**2)
    constrained by given observed directional moments,

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array with shape [number_of_frequencies,number_of_direction]
    representing the directional distribution of the waves at each frequency.

    Minimize the variance of the solution:

           integrate D**2 over directions

    such that the resulting distribution D reproduces the observed moments.

    Implementation notes:
    - we formulate the problem as a standard Quadratic Programming problem which
      can them be solved efficiently with the qpsolvers package.
    """

    a1 = numpy.atleast_1d(a1)
    b1 = numpy.atleast_1d(b1)
    a2 = numpy.atleast_1d(a2)
    b2 = numpy.atleast_1d(b2)

    number_of_frequencies = len(a1)
    directional_distribution = numpy.zeros(
        (number_of_frequencies, len(directions_radians))
    )

    constraint_matrix = get_constraint_matrix(directions_radians)
    rhs = get_rhs(a1, b1, a2, b2)
    identity_matrix = numpy.diag(numpy.ones_like(directions_radians), 0)

    zeros = numpy.zeros_like(directions_radians)
    direction_increment = get_direction_increment(directions_radians)
    upperbound = numpy.ones_like(directions_radians) / direction_increment.min()

    for ifreq in range(0, number_of_frequencies):
        res = solve_ls(
            R=identity_matrix,
            s=zeros,
            # minimizing |Rx-b|**2
            lb=zeros,
            ub=upperbound,
            # lb: non-negative; ub: binwidth * ub = 1
            A=constraint_matrix,
            b=rhs[ifreq, :],
            verbose=False
            # with hard constraint that Ax=b
        )

        if res is None:
            raise Exception("No solution")

        directional_distribution[ifreq, :] = res

    directional_distribution[directional_distribution < 0] = 0
    return numpy.squeeze(directional_distribution)


### mem.py
import numpy
from numba import njit


def mem(
    directions_radians: numpy.ndarray,
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
    progress,
    **kwargs
) -> numpy.ndarray:

    number_of_frequencies = a1.shape[-1]
    number_of_points = a1.shape[0]

    directional_distribution = numpy.zeros(
        (number_of_points, number_of_frequencies, len(directions_radians))
    )

    for ipoint in range(0, number_of_points):
        progress.update(1)
        directional_distribution[ipoint, :, :] = _mem(
            directions_radians,
            a1[ipoint, :],
            b1[ipoint, :],
            a2[ipoint, :],
            b2[ipoint, :],
        )
    return directional_distribution


def _mem(
    directions_radians: numpy.ndarray,
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
) -> numpy.ndarray:
    """
    This function uses the maximum entropy method by Lygre and Krogstadt (1986,JPO)
    to estimate the directional shape of the spectrum. Enthropy is defined in the
    Boltzmann sense (log D)

    Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the directional
    distribution in ocean wave spectra. Journal of Physical Oceanography, 16(12), 2052-2060.

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]. (going to, anti-clockswise from east)

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array with shape [number_of_frequencies,number_of_direction]
    representing the directional distribution of the waves at each frequency.

    Maximize the enthrophy of the solution with entrophy defined as:

           integrate log(D) over directions

    such that the resulting distribution D reproduces the observed moments.

    :return: Directional distribution as a numpy array

    Note that:
    d1 = a1; d2 =b1; d3 = a2 and d4=b2 in the defining equations 10.
    """

    number_of_directions = len(directions_radians)

    c1 = a1 + 1j * b1
    c2 = a2 + 1j * b2
    #
    # Eq. 13 L&K86
    #
    Phi1 = (c1 - c2 * numpy.conj(c1)) / (1 - c1 * numpy.conj(c1))
    Phi2 = c2 - Phi1 * c1
    #
    e1 = numpy.exp(-directions_radians * 1j)
    e2 = numpy.exp(-directions_radians * 2j)

    numerator = 1 - Phi1 * numpy.conj(c1) - Phi2 * numpy.conj(c2)
    denominator = (
        numpy.abs(1 - Phi1[:, None] * e1[None, :] - Phi2[:, None] * e2[None, :]) ** 2
    )

    D = numpy.real(numerator[:, None] / denominator) / numpy.pi / 2

    # Normalize to 1. in discrete sense
    integralApprox = numpy.sum(D, axis=-1) * numpy.pi * 2.0 / number_of_directions
    D = D / integralApprox[:, None]

    return numpy.squeeze(D)


@njit(cache=True)
def numba_mem(
    directions_radians: numpy.ndarray,
    a1: float,
    b1: float,
    a2: float,
    b2: float,
) -> numpy.ndarray:
    """
    This function uses the maximum entropy method by Lygre and Krogstadt (1986,JPO)
    to estimate the directional shape of the spectrum. Enthropy is defined in the
    Boltzmann sense (log D)

    Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the directional
    distribution in ocean wave spectra. Journal of Physical Oceanography, 16(12), 2052-2060.

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]. (going to, anti-clockswise from east)

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array with shape [number_of_frequencies,number_of_direction]
    representing the directional distribution of the waves at each frequency.

    Maximize the enthrophy of the solution with entrophy defined as:

           integrate log(D) over directions

    such that the resulting distribution D reproduces the observed moments.

    :return: Directional distribution as a numpy array

    Note that:
    d1 = a1; d2 =b1; d3 = a2 and d4=b2 in the defining equations 10.
    """

    number_of_directions = len(directions_radians)

    c1 = a1 + 1j * b1
    c2 = a2 + 1j * b2
    #
    # Eq. 13 L&K86
    #
    Phi1 = (c1 - c2 * numpy.conj(c1)) / (1 - c1 * numpy.conj(c1))
    Phi2 = c2 - Phi1 * c1
    #
    e1 = numpy.exp(-directions_radians * 1j)
    e2 = numpy.exp(-directions_radians * 2j)

    numerator = 1 - Phi1 * numpy.conj(c1) - Phi2 * numpy.conj(c2)
    denominator = numpy.abs(1 - Phi1 * e1 - Phi2 * e2) ** 2

    D = numpy.real(numerator / denominator) / numpy.pi / 2

    # Normalize to 1. in discrete sense
    integralApprox = numpy.sum(D, axis=-1) * numpy.pi * 2.0 / number_of_directions
    D = D / integralApprox

    return D



### mem2.py 
"""
Implementation of the "MEM2" method:

see Kim1995:

    Kim, T., Lin, L. H., & Wang, H. (1995). Application of maximum entropy method
    to the real sea data. In Coastal Engineering 1994 (pp. 340-355).

    link: https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/4967/4647
    (working as of May 29, 2022)

and references therein.

"""
import numpy
from scipy.optimize import root
import typing
#from roguewave.wavespectra.estimators.utils import get_direction_increment
#from roguewave.wavespectra.estimators.mem import numba_mem
from numba import njit, prange
from numba.typed import Dict as NumbaDict
from numba.core import types
from numpy.linalg import norm
from numba_progress import ProgressBar

# Settings for numba JIT compilation- whether to use fast math and parallel optimizations when possible.
_FASTMATH = True
_PARALLEL = False

# Numerical settings used in solving for the mem2 distribution
NUMERICS = {
    # absolute tolerence stopping criterium. let moment = [ a1,b1,a2,b2] and let iterate_moment contain the moments
    # calculated from the current estmitaed distribution. The stopping criterium is:
    #     norm( moment-iterate_moment ) < atol
    "atol": 0.01,
    # Maximum number of iterations
    "max_iter": 100,
    # Maximum number of subiterations in the line search algorithm. Typically deep line search activates only when
    # the convergence is poor anyway.
    "max_line_search_depth": 8,
    # If we fall back to least squares estimate of the newton update we have an ill-conditioned system, and solve
    # the system approximately removing the smallest singular values. rcond it the ration of smallest divided by largest
    # singular value.
    "rcond": 1e-6,
    # Convergence is mostly (based on limited testing) poor for narrow distributions (large lagrange multipliers). If
    # we fail to converge we fall back to the mem estimate which has no such issues. For narrow distributions this is
    # hopefully fine.
    "use_mem_when_failing_to_converge": True,
}

# Entry Function
# =============================================================================


def mem2(
    directions_radians: numpy.ndarray,
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
    progress_bar: ProgressBar = None,
    solution_method="newton",
    solver_config=None,
) -> numpy.ndarray:
    """

    :param directions_radians:
    :param a1:
    :param b1:
    :param a2:
    :param b2:
    :param solution_method:
    :return:
    """

    if solver_config is None:
        solver_config = NUMERICS

    else:
        solver_config = NUMERICS | solver_config

    if solution_method == "scipy":
        func = mem2_scipy_root_finder
        kwargs = {}

    elif solution_method == "newton":
        func = mem2_newton
        numba_solver_config = NumbaDict.empty(
            key_type=types.unicode_type, value_type=types.float64
        )
        for key in solver_config:
            numba_solver_config[key] = solver_config[key]

        kwargs = {"config": numba_solver_config}

    elif solution_method == "approximate":
        func = mem2_newton
        kwargs = {"approximate": True}

    else:
        raise ValueError("Unknown method")

    return func(directions_radians, a1, b1, a2, b2, progress_bar, **kwargs)


# Scipy Implementation
# =============================================================================
def mem2_scipy_root_finder(
    directions_radians: numpy.ndarray,
    a1: typing.Union[numpy.ndarray, float],
    b1: typing.Union[numpy.ndarray, float],
    a2: typing.Union[numpy.ndarray, float],
    b2: typing.Union[numpy.ndarray, float],
    progress,
    **kwargs
) -> numpy.ndarray:
    """
    Return the directional distribution that maximizes Shannon [ - D log(D) ]
    enthrophy constrained by given observed directional moments,

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array with shape [number_of_frequencies,number_of_direction]
    representing the directional distribution of the waves at each frequency.

    Maximize the enthrophy of the solution with entrophy defined as:

           integrate - D * log(D) over directions

    such that the resulting distribution D reproduces the observed moments.

    """

    number_of_frequencies = a1.shape[-1]
    number_of_points = a1.shape[0]

    directional_distribution = numpy.zeros(
        (number_of_points, number_of_frequencies, len(directions_radians))
    )

    direction_increment = get_direction_increment(directions_radians)

    twiddle_factors = numpy.empty((4, len(directions_radians)))
    twiddle_factors[0, :] = numpy.cos(directions_radians)
    twiddle_factors[1, :] = numpy.sin(directions_radians)
    twiddle_factors[2, :] = numpy.cos(2 * directions_radians)
    twiddle_factors[3, :] = numpy.sin(2 * directions_radians)

    guess = initial_value(a1, b1, a2, b2)
    for ipoint in range(0, number_of_points):
        progress.update(1)
        for ifreq in range(0, number_of_frequencies):
            #
            moments = numpy.array(
                [
                    a1[ipoint, ifreq],
                    b1[ipoint, ifreq],
                    a2[ipoint, ifreq],
                    b2[ipoint, ifreq],
                ]
            )

            if numpy.any(numpy.isnan(guess[ipoint, ifreq, :])):
                continue

            res = root(
                moment_constraints,
                guess[ipoint, ifreq, :],
                args=(twiddle_factors, moments, direction_increment),
                method="lm",
            )
            lambas = res.x

            directional_distribution[ipoint, ifreq, :] = mem2_directional_distribution(
                lambas, direction_increment, twiddle_factors
            )

    return directional_distribution


# Numba Implementation
# =============================================================================


# spatial iteration
# ---------------------


# To note; enabling caching seems to not play nice with paralel
@njit(parallel=_PARALLEL, cache=(not _PARALLEL))
def mem2_newton(
    directions_radians: numpy.ndarray,
    a1: numpy.ndarray,
    b1: numpy.ndarray,
    a2: numpy.ndarray,
    b2: numpy.ndarray,
    progress_bar: ProgressBar = None,
    config: NumbaDict = None,
    approximate: bool = False,
) -> numpy.ndarray:
    """
    Return the directional distribution that maximizes Shannon [ - D log(D) ]
    enthrophy constrained by given observed directional moments.

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]

    :param a1: 1d array of cosine directional moment as function of position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param b1: 1d array of sine directional moment as function of position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param a2: 1d array of double angle cosine directional moment as function of position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param b2: 1d array of double angle sine directional moment as function of position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param progress_bar: Progress bar instance if updates are desired.

    :return: array with shape [numbrt_of_points, number_of_frequencies,number_of_direction]
    representing the directional distribution of the waves at each frequency.

    Maximize the enthrophy of the solution with entrophy defined as:

           integrate - D * log(D) over directions

    such that the resulting distribution D reproduces the observed moments.

    """

    number_of_frequencies = a1.shape[-1]
    number_of_points = a1.shape[0]

    directional_distribution = numpy.zeros(
        (number_of_points, number_of_frequencies, len(directions_radians))
    )

    direction_increment_downward_difference = (
        directions_radians - numpy.roll(directions_radians, 1) + numpy.pi
    ) % (2 * numpy.pi) - numpy.pi

    direction_increment_upward_difference = (
        -(directions_radians - numpy.roll(directions_radians, -1) + numpy.pi)
        % (2 * numpy.pi)
        - numpy.pi
    )

    direction_increment = (
        direction_increment_downward_difference + direction_increment_upward_difference
    ) / 2

    # Calculate the needed Fourier transform twiddle factors to calculate moments.
    twiddle_factors = numpy.empty((4, len(directions_radians)))
    twiddle_factors[0, :] = numpy.cos(directions_radians)
    twiddle_factors[1, :] = numpy.sin(directions_radians)
    twiddle_factors[2, :] = numpy.cos(2 * directions_radians)
    twiddle_factors[3, :] = numpy.sin(2 * directions_radians)

    guess = initial_value(a1, b1, a2, b2)
    for ipoint in prange(0, number_of_points):
        if progress_bar is not None:
            progress_bar.update(1)

        # Note; entries to directional_distribution[ipoint, :, :] is modified in the call below. This avoids creation
        # of memory for the resulting array at the expense of allowing for side-effects.
        _mem2_newton_point(
            directional_distribution[ipoint, :, :],
            a1[ipoint, :],
            b1[ipoint, :],
            a2[ipoint, :],
            b2[ipoint, :],
            guess[ipoint, :, :],
            direction_increment,
            twiddle_factors,
            config,
            approximate,
        )

    return directional_distribution


# frequency iteration
# ----------------------


@njit(cache=True)
def _mem2_newton_point(
    out,
    a1,
    b1,
    a2,
    b2,
    guess,
    direction_increment,
    twiddle_factors,
    config=None,
    approximate=False,
):
    """

    :param out: a (view) of the array that will containt the output
    :param a1: 1d array of cosine directional moment as function of frequency,
    :param b1: 1d array of sine directional moment as function of frequency,
    :param a2: 1d array of double angle cosine directional moment as function of frequency,
    :param b2: 1d array of double angle sine directional moment as function of frequency,
    :param guess: initial guess of the lagrange multipliers
    :param direction_increment: directional stepsize used in the integration, nd-array
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
    :param config: numerical settings, see description at NUMERICS at top of file.
    :param approximate: whether or not to use the approximate relations.
    :return: None - we use side-effects to pass the results back to the caller (modifying out)
    """
    number_of_frequencies = a1.shape[0]
    for ifreq in range(0, number_of_frequencies):
        #
        moments = numpy.array([a1[ifreq], b1[ifreq], a2[ifreq], b2[ifreq]])
        out[ifreq, :] = mem2_newton_solver(
            moments,
            guess[ifreq, :],
            direction_increment,
            twiddle_factors,
            config,
            approximate,
        )


# mem2 numerical solver
# ----------------------


@njit(cache=True)
def mem2_newton_solver(
    moments: numpy.ndarray,
    guess: numpy.ndarray,
    direction_increment: numpy.ndarray,
    twiddle_factors: numpy.ndarray,
    config=None,
    approximate=False,
) -> numpy.ndarray:
    """
    Newton iteration to find the solution to the non-linear system of constraint equations defining the lagrange
    multipliers in the MEM2 method. Because the Lagrange multipliers enter the equations as exponents the system can
    be unstable to solve numerically.

    :param moments: the normalized directional moments [a1,b1,a2,b2]
    :param guess: first guess for the lagrange multipliers (ndarray, length 4)
    :param direction_increment: directional stepsize used in the integration, nd-array
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
    :param config: numerical settings, see description at NUMERICS at top of file.
    :param approximate: whether or not to use the approximate relations.
    :return:
    """
    if config is None:
        max_iter = 100
        rcond = 1e-6
        atol = 0.01
        max_line_search_depth = 8
        use_mem_when_failing_to_converge = True

    else:
        max_iter = config["max_iter"]
        rcond = config["rcond"]
        atol = config["atol"]
        max_line_search_depth = config["max_line_search_depth"]
        use_mem_when_failing_to_converge = (
            config["use_mem_when_failing_to_converge"] > 0.0
        )

    directional_distribution = numpy.empty(len(direction_increment))
    if numpy.any(numpy.isnan(guess)):
        directional_distribution[:] = 0
        return directional_distribution

    if approximate:
        directional_distribution[:] = mem2_directional_distribution(
            guess, direction_increment, twiddle_factors
        )
        return directional_distribution

    current_iterate = guess
    current_func = moment_constraints(
        current_iterate,
        twiddle_factors,
        moments,
        direction_increment,
    )

    jacobian = numpy.empty((4, 4))

    convergence = False
    for iter in range(0, max_iter):

        # Stopping criterium
        magnitude_cur_func_eval = norm(current_func)
        if magnitude_cur_func_eval < atol:
            convergence = True
            break

        #
        # Compute jacobian, and find newton iterate innovation as we solve for:
        #
        #       jacobian @ delta = - current_iterate_func_eval
        #
        # with:
        #
        #       delta = next_lagrange_multiplier_iterate-cur_lagrange_multiplier_iterate

        jacobian = mem2_jacobian(
            current_iterate, twiddle_factors, direction_increment, jacobian
        )
        try:
            update_iterate = solve_cholesky(jacobian, -current_func)
        except Exception:
            update_iterate = numpy.linalg.lstsq(jacobian, -current_func, rcond=rcond)[0]

        magnitude_current_iterate = norm(current_iterate)
        magnitude_update = norm(update_iterate)

        # Do a line search for the optimum decrease. This is intended to stabilize the algorithm
        # as the equations are ill-posed.
        line_search_factor = 1
        for ii in range(max_line_search_depth):
            next_iterate = current_iterate + line_search_factor * update_iterate
            next_func = moment_constraints(
                next_iterate, twiddle_factors, moments, direction_increment
            )

            if norm(next_func) < magnitude_cur_func_eval:
                # If we are decreasing- continue
                current_func = next_func
                current_iterate = next_iterate
                break
            else:
                # The update may be too big as we are not decreasing the cost function magnitude. We will decrease the
                # step size we take - but keep the direction of the step the same.
                inverse_relative_update = magnitude_current_iterate / magnitude_update
                line_search_factor = min(
                    inverse_relative_update, line_search_factor / 2
                )
        else:
            # The linesearch failed. We could not find a factor that ensures the next function estimate is closer
            # to 0.
            convergence = False
            break
    else:
        # We failed to converge after the maximum number of iterations.
        convergence = False

    if not convergence:
        if use_mem_when_failing_to_converge:
            directions = numpy.arctan2(twiddle_factors[1, :], twiddle_factors[0, :])
            directional_distribution[:] = numba_mem(
                directions, moments[0], moments[1], moments[2], moments[3]
            )
        else:
            raise ValueError("we did not converge")

    directional_distribution[:] = mem2_directional_distribution(
        current_iterate, direction_increment, twiddle_factors
    )

    return directional_distribution


# mem2 functions
# ----------------------


@njit(cache=True, fastmath=_FASTMATH)
def moment_constraints(lambdas, twiddle_factors, moments, direction_increment):
    """
    Construct the nonlinear equations we need to solve for lambda. The constrainst are the difference between the
    desired moments a1,b1,a2,b2 and the moment calculated from the current distribution guess and for a perfect fit
    should be 0.

    To note: we differ from Kim et al here who formulate the constraints using unnormalized equations. Here we opt to
    use the normalized version as that allows us to cast the error / or mismatch directly in terms of an error in the
    moments.

    :param lambdas: the lagrange multipliers
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
    :param moments: [a1,b1,a2,b2]
    :param direction_increment: directional stepsize used in the integration, nd-array
    :return: array (length=4) with the difference between desired moments and those calculated from the current
        approximate distribution
    """

    # Get the current estimate of the directional distribution
    dist = mem2_directional_distribution(lambdas, direction_increment, twiddle_factors)
    out = numpy.zeros(4)
    for mm in range(0, 4):
        # note - the part after the "-" is just a discrete approximation of the Fourier sine/cosine amplitude (moment)
        out[mm] = moments[mm] - numpy.sum(
            (twiddle_factors[mm, :]) * dist * direction_increment
        )

    return out


@njit(cache=True, fastmath=_FASTMATH)
def mem2_jacobian(lagrange_multiplier, twiddle_factors, direction_increment, jacobian):
    """
    Calculate the jacobian of the constraint equations. The resulting jacobian is a square and positive definite matrix

    :param lambdas: the lagrange multipliers
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
    :param direction_increment: directional stepsize used in the integration, nd-array

    :return: a 4 by 4 matrix that is the Jacobian of the constraint equations.
    """
    inner_product = numpy.zeros(twiddle_factors.shape[1])
    for jj in range(0, 4):
        inner_product = inner_product + lagrange_multiplier[jj] * twiddle_factors[jj, :]

    # We subtract the minimum to ensure that the values in the exponent do not become too large. This amounts to
    # multiplyig with a constant - which is fine since we normalize anyway. Effectively- this avoids overflow errors
    # (or infinities) - at the expense of underflowing (which is less of an issue).
    #
    inner_product = inner_product - numpy.min(inner_product)

    normalization = 1 / numpy.sum(numpy.exp(-inner_product) * direction_increment)
    shape = numpy.exp(-inner_product)

    normalization_derivative = numpy.zeros(4)
    for mm in range(0, 4):
        normalization_derivative[mm] = normalization * numpy.sum(
            twiddle_factors[mm, :] * numpy.exp(-inner_product) * direction_increment
        )

    # To note- we have to multiply seperately to avoid potential underflow/overflow errors.
    normalization_derivative = normalization_derivative * normalization

    shape_derivative = numpy.zeros((4, twiddle_factors.shape[1]))
    for mm in range(0, 4):
        shape_derivative[mm, :] = -twiddle_factors[mm, :] * shape

    for mm in range(0, 4):
        # we make use of symmetry and only explicitly calculate up to the diagonal
        for nn in range(0, mm + 1):
            jacobian[mm, nn] = -numpy.sum(
                twiddle_factors[mm, :]
                * direction_increment
                * (
                    normalization * shape_derivative[nn, :]
                    + shape * normalization_derivative[nn]
                ),
                -1,
            )
            if nn != mm:
                jacobian[nn, mm] = jacobian[mm, nn]
    return jacobian


@njit(cache=True, fastmath=_FASTMATH)
def mem2_directional_distribution(
    lagrange_multiplier, direction_increment, twiddle_factors
) -> numpy.ndarray:
    """
    Given the solution for the Lagrange multipliers- reconstruct the directional
    distribution.
    :param lagrange_multiplier: the lagrange multipliers
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
    :param direction_increment: directional stepsize used in the integration, nd-array
    :return: Directional distribution arrasy as a function of directions
    """
    inner_product = numpy.zeros(twiddle_factors.shape[1])
    for jj in range(0, 4):
        inner_product = inner_product + lagrange_multiplier[jj] * twiddle_factors[jj, :]

    inner_product = inner_product - numpy.min(inner_product)

    normalization = 1 / numpy.sum(numpy.exp(-inner_product) * direction_increment)
    return numpy.exp(-inner_product) * normalization


@njit(cache=True, fastmath=_FASTMATH)
def initial_value(
    a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray, b2: numpy.ndarray
):
    """
    Initial guess of the Lagrange Multipliers according to the "MEM AP2" approximation
    found im Kim1995

    :param a1: moment a1
    :param b1: moment b1
    :param a2: moment a2
    :param b2: moment b2
    :return: initial guess of the lagrange multipliers, with the same leading dimensions as input.
    """
    guess = numpy.empty((*a1.shape, 4))
    fac = 1 + a1**2 + b1**2 + a2**2 + b2**2
    guess[..., 0] = 2 * a1 * a2 + 2 * b1 * b2 - 2 * a1 * fac
    guess[..., 1] = 2 * a1 * b2 - 2 * b1 * a2 - 2 * b1 * fac
    guess[..., 2] = a1**2 - b1**2 - 2 * a2 * fac
    guess[..., 3] = 2 * a1 * b1 - 2 * b2 * fac
    return guess


@njit(cache=True, fastmath=True)
def solve_cholesky(matrix, rhs):
    """
    Solve using cholesky decomposition according to the Choleskyâ€“Banachiewicz algorithm.
    See: https://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky_algorithm
    """
    M, N = matrix.shape
    x = numpy.zeros(M)
    cholesky_decomposition = numpy.zeros((M, M))
    inv = numpy.zeros(M)

    for mm in range(0, M):
        forward_sub_sum = rhs[mm]
        for nn in range(0, mm):
            sum = matrix[mm, nn]
            for kk in range(0, nn):
                sum -= cholesky_decomposition[mm, kk] * cholesky_decomposition[nn, kk]

            cholesky_decomposition[mm, nn] = inv[nn] * sum
            forward_sub_sum += -cholesky_decomposition[mm, nn] * x[nn]

        sum = matrix[mm, mm]
        for kk in range(0, mm):
            sum -= cholesky_decomposition[mm, kk] ** 2

        if sum <= 0.0:
            raise ValueError(
                "Matrix not positive definite, likely due to finite precision errors."
            )

        cholesky_decomposition[mm, mm] = numpy.sqrt(sum)
        inv[mm] = 1 / cholesky_decomposition[mm, mm]
        x[mm] = forward_sub_sum * inv[mm]

    # Backward Substitution (in place)
    for mm in range(0, M):
        kk = M - mm - 1
        sum = x[kk]
        for nn in range(kk + 1, N):
            sum += -cholesky_decomposition[nn, kk] * x[nn]
        x[kk] = sum * inv[kk]
    return x


### utils.py
import numpy

def get_direction_increment(
        directions_radians: numpy.ndarray) -> numpy.ndarray:
    """
    calculate the stepsize used for midpoint integration. The directions
    represent the center of the interval - and we want to find the dimensions of
    the interval (difference between the preceeding and succsesive midpoint).

    :param directions_radians: array of radian directions
    :return: array of radian intervals
    """

    # Calculate the forward difference appending the first entry to the back
    # of the array. Use modular trickery to ensure the angle is in [-pi,pi]
    forward_diff = (numpy.diff(directions_radians,
                               append=directions_radians[0]) + numpy.pi) % (
                           2 * numpy.pi) - numpy.pi

    # Calculate the backward difference prepending the last entry to the front
    # of the array. Use modular trickery to ensure the angle is in [-pi,pi]
    backward_diff = (numpy.diff(directions_radians,
                                prepend=directions_radians[-1]) + numpy.pi) % (
                            2 * numpy.pi) - numpy.pi

    # The interval we are interested in is the average of the forward and backward
    # differences.
    return (forward_diff + backward_diff) / 2


def get_constraint_matrix(directions_radians: numpy.ndarray) -> numpy.ndarray:
    """
    Define the matrix M that can be used in the matrix product M@D (with D the
    directional distribution) such that:

            M@D = [1,a1,b1,a2,b2]^T

    with a1,b1 etc the directional moments at a given frequency.

    :param directions_radians: array of radian directions
    :return:
    """
    number_of_dir = len(directions_radians)
    constraints = numpy.zeros((5, number_of_dir))
    direction_increment = get_direction_increment(directions_radians)
    constraints[0, :] = direction_increment
    constraints[1, :] = direction_increment * numpy.cos(directions_radians)
    constraints[2, :] = direction_increment * numpy.sin(directions_radians)
    constraints[3, :] = direction_increment * numpy.cos(2 * directions_radians)
    constraints[4, :] = direction_increment * numpy.sin(2 * directions_radians)
    return constraints


def get_rhs(a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray,
            b2: numpy.ndarray) -> numpy.ndarray:
    """
    Define the matrix rhs that for each row contains the directional moments
    at a given frequency:

    rhs = [ 1, a1[0],b1[0],a2[0],b2[0],
            |    |    |      |    |
            N, a1[0],b1[0],a2[0],b2[0] ]

    These rows are use as the "right hand side" in the linear constraints
    (see get_constraint_matrix)

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array ( number of frequencies by 5) that for each row contains
    the directional moments at a given frequency
    """
    rhs = numpy.array([numpy.ones_like(a1), a1, b1, a2, b2]).transpose()
    return rhs
