"""
precip
======

computation of bulk precipitation characteristics

.. autosummary::
    :toctree: generated/

    compute_lwc
    compute_elwc
    compute_rainfall_rate
    compute_breakup_prob
    compute_dsd_shed_water
    compute_dsd_breakup_water
    psd_ic_func
    psd_hail

"""

from copy import deepcopy

import numpy as np

from pytmatrix.psd import UnnormalizedGammaPSD, ExponentialPSD

from part_descrip import compute_axis_ratio_ic


def compute_lwc(diam, delta_d, psd_vals):
    """
    Computes the liquid water content

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    delta_d : array of floats
        diameter bin size (mm)
    psd_vals : array of floats
        DSD value at each diameter bin (number of particles)

    Returns
    -------
    LWC : array of floats
       the liquid water content (mm3/m3)

    """
    return np.pi/6.*np.sum(np.power(diam, 3.)*psd_vals*delta_d)


def compute_elwc(delta_d, mass, psd_vals, dens_w=0.0009999720):
    """
    Computes the equivalent liquid water content

    Parameters
    ----------
    delta_d : array of floats
        diameter bin size (mm)
    mass : array of floats
        particle mass (g)
    psd_vals : array of floats
        PSD value at each diameter bin (number of particles)
    dens_w : float
        water density at 4 deg C (g/mm3)

    Returns
    -------
    ELWC : array of floats
       the liquid water content (mm3/m3)

    """

    return np.sum(mass*psd_vals*delta_d)/dens_w


def compute_rainfall_rate(diam, delta_d, psd_vals, vel):
    """
    Computes the rainfall rate

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    delta_d : array of floats
        diameter bin size (mm)
    psd_vals : array of floats
        DSD value at each diameter bin (number of particles)
    vel : array of floats
        terminal velocity of the raindrops

    Returns
    -------
    RR : array of floats
       the rainfall rate (mm/h)

    """
    return 0.6*np.pi/6.*1e-3*np.sum(vel*np.power(diam, 3.)*psd_vals*delta_d)


def compute_equi_rainfall_rate(delta_d, mass, psd_vals, vel,
                               dens_w=0.0009999720):
    """
    Computes the equivalent rainfall rate

    Parameters
    ----------
    delta_d : array of floats
        diameter bin size (mm)
    mass : array of floats
        hydrometeor mass (g)
    psd_vals : array of floats
        PSD value at each diameter bin (number of particles)
    vel : array of floats
        terminal velocity of the hydrometeor (m/s)
    dens_w : float
        water density at 4 deg C (g/mm3)

    Returns
    -------
    eRR : array of floats
       the equivalent rainfall rate (mm/h)

    """

    return np.sum(mass*vel*psd_vals*delta_d)/dens_w*3.6/1000.


def compute_breakup_prob(h0m, prob_break, fmw, d_ext, alt, h_b=400.,
                         d_breakup=7.95):
    """
    Computes the probability that a large drop breaks up

    Parameters
    ----------
    h0m : array of floats
        height at which a hailstone of a given initial size melts and becomes
        a raindrop 8 mm or larger (masl)
    prob_break : array of floats
        the probability of a breakup at altitude right above the current
        altitude
    fmw : array of floats
        fraction mass water
    d_ext : array of floats
        external hailstone diameter
    alt : float
        altitude at which the probability of break up is evaluated (masl)
    h_b : float
        altitude parameter
    d_breakup : float
        diameter of drops that can breakup (mm)

    Returns
    -------
    h0m_out : array of floats
       the updated h0m
    prob_break_out : array of floats
        the updated breakup probability

    """
    h0m_out = deepcopy(h0m)
    prob_break_out = deepcopy(prob_break)

    ind_break = np.where(
        (np.isclose(fmw, 1.)) & (d_ext >= d_breakup) & (np.isnan(h0m_out)))[0]
    if ind_break.size > 0:
        h0m_out[ind_break] = alt
    prob_break_out[~np.isnan(h0m_out)] = 1.-np.exp(
        -np.power((h0m_out[~np.isnan(h0m_out)]-alt)/(1.2*h_b), 2.))

    return h0m_out, prob_break_out


def compute_dsd_shed_water(diam_h, psd_vals_hail, diam_rd, mass_shed,
                           dens_w=0.0009999720, lamb=2., mu=2.):
    """
    Computes the DSD of water shedded in the process of melting a hailstone

    Parameters
    ----------
    diam_h : array of floats
        equivalent volume diameter of the mother hailstone (mm)
    psd_vals_hail : array of floats
        number of hailstones for each size
    diam_rd : array of floats
        equivalent volume diameter of the raindrops formed by shed water (mm)
    mass_shed : array of floats
        mass of water shed by each hailstone size (g)
    dens_w : float
        water density at 4 deg C (g/mm3)
    lamb, mu : float
        parameters of the gamma distribution DSD of the shed water

    Returns
    -------
    psd_func : function
       Function describing the DSD of shed water

    """
    # total mass of shed water
    delta_d_hail = diam_h-np.append(0, diam_h[:-1])
    mass_shed_t = np.sum(mass_shed*psd_vals_hail*delta_d_hail)

    delta_d_drop = diam_rd-np.append(0, diam_rd[:-1])
    mass_rd = dens_w*np.pi/6.*np.power(diam_rd, 3.)
    n0 = mass_shed_t/np.sum(
        mass_rd*np.power(diam_rd, mu)*np.exp(-lamb*diam_rd)*delta_d_drop)

    return UnnormalizedGammaPSD(N0=n0, Lambda=lamb, mu=mu, D_max=diam_rd[-1])


def compute_dsd_breakup_water(diam_h, psd_breakup_hail, diam_rd, water_mass,
                              dens_w=0.0009999720, lamb=0.453):
    """
    Computes the DSD of drops resulting from break up of large drops created
    in the process of melting a hailstone

    Parameters
    ----------
    diam_h : array of floats
        equivalent volume diameter of the mother hailstone (mm)
    psd_vals_hail : array of floats
        number of hailstones for each size bin
    diam_rd : array of floats
        equivalent volume diameter of the raindrops formed by breakup (mm)
    water_mass : array of floats
        mass of totally melted hailstones
    dens_w : float
        water density at 4 deg C (g/mm3)
    lamb : float
        parameters of the exponential distribution of the breakup drops

    Returns
    -------
    psd_func : function
       Function describing the DSD of breakup water

    """
    # total mass of breakup water
    delta_d_hail = diam_h-np.append(0, diam_h[:-1])
    mass_break_t = np.sum(water_mass*psd_breakup_hail*delta_d_hail)

    delta_d_drop = diam_rd-np.append(0, diam_rd[:-1])
    mass_rd = dens_w*np.pi/6.*np.power(diam_rd, 3.)
    n0 = mass_break_t/np.sum(mass_rd*np.exp(-lamb*diam_rd)*delta_d_drop)

    return ExponentialPSD(N0=n0, Lambda=lamb)


def psd_ic_func(length, ar_coeff, ar_exponent, n0, lamb):
    """
    Computes the PSD of ice crystals as a function of maximum dimension

    Parameters
    ----------
    length : float
        maximum dimension (mm)
    ar_coeff, ar_exponent : float
        coefficient and exponent of the power law relating axis ratio and
        maximum dimension
    N0, lamb : float
        parameters of the exponential distribution

    Returns
    -------
    psd_val : float
       PSD value

    """
    axis_ratio = compute_axis_ratio_ic(
        np.array([length]), ar_coeff, ar_exponent)[0]
    diam = length*np.power(axis_ratio, 1./3.)
    psd_func = ExponentialPSD(N0=n0, Lambda=lamb)

    return psd_func(diam)


def psd_hail(diam, n_g=8000., lamb_g=1.6, d_max_g=8., d_max_h=35,
             lamb_h=0.27, coeff_nh=800.):
    """
    Computes hail PSD. Hail PSD is a bi-modal exponential

    Parameters
    ----------
    diam : array of float
        hailstone diameter
    n_g, lamb_g, d_max_g: float
        Parameters of the graupel part of the PSD
    d_max_h, lamb_h, coeff_nh : float
        Parameters to compute the hail part of the PSD

    Returns
    -------
    psd_h_vals: array of floats
        the values of the psd for each diam

    """
    # PSD graupel
    psd_g = ExponentialPSD(N0=n_g, Lambda=lamb_g, D_max=d_max_g)

    # PSD hail
    n_h = coeff_nh*np.power(lamb_h, 4.11)
    psd_h = ExponentialPSD(N0=n_h, Lambda=lamb_h, D_max=d_max_h)
    return psd_h(diam) + psd_g(diam)
