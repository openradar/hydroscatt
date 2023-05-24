"""
part_descrip
============

functions to compute parameters describing particle properties

.. autosummary::
    :toctree: generated/

    compute_axis_ratio_brandes
    compute_axis_ratio_mh
    compute_axis_ratio_dh
    compute_axis_ratio_snow
    compute_axis_ratio_melting_snow
    compute_axis_ratio_ic
    compute_equi_vol_diam
    compute_ms_core_diameter
    compute_volume
    compute_ice_mass
    compute_snow_mass
    compute_hail_water_max
    compute_ice_crystal_density
    compute_snow_density
    compute_melting_snow_density
    compute_hailstone_density
    compute_hs_air_vol
    compute_velocity_rain
    compute_velocity_ic
    compute_velocity_snow
    compute_velocity_melting_snow
    compute_velocity_melting_hail
    compute_reynolds_number
    compute_reynolds_number_snow
    compute_snow_melting
    compute_hail_melting
    snowflake_change
    hailstone_change

"""

# from warnings import warn
from copy import deepcopy

import numpy as np
# g : gravitational acceleration 9.80665 m/s2
# gas_constant: 8.314462618 J/(mol K)
from scipy.constants import g, gas_constant

from pytmatrix.tmatrix_aux import dsr_thurai_2007


def compute_axis_ratio_brandes(diam):
    """
    Computes the axis ratio of raindrops using the relation of Brandes et al
    (2002)

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)

    Returns
    -------
    axis_ratio : array of floats
       axis ratio of raindrops

    """
    axis_ratio = (
        0.9951+0.02510*diam-0.03544*diam*diam
        + 0.005303*np.power(diam, 3.)-0.0002492*np.power(diam, 4.))
    axis_ratio[axis_ratio > 1.] = 1.

    return axis_ratio


def compute_axis_ratio_mh(diam, fmw):
    """
    Computes the axis ratio of melting hailstones

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    fmw : array of floats
        fraction mass of water

    Returns
    -------
    axis_ratio_mh : array of floats
       axis ratio of a melting hailstone

    """
    axis_ratio_mh = np.zeros(diam.size)
    axis_ratio_w = np.zeros(diam.size)

    for ind, diam_part in enumerate(diam):
        axis_ratio_w[ind] = dsr_thurai_2007(diam_part)
    # limit the axis ratio of big drops
    axis_ratio_w[axis_ratio_w < 0.56] = 0.56

    axis_ratio_h = compute_axis_ratio_dh(diam)

    ind = np.where(fmw < 0.2)[0]
    if ind.size > 0:
        axis_ratio_mh[ind] = (
            axis_ratio_h[ind]-5.*(axis_ratio_h[ind]-0.8)*fmw[ind])
    ind = np.where((fmw >= 0.2) & (fmw <= 0.8))[0]
    if ind.size > 0:
        axis_ratio_mh[ind] = 0.88-0.4*fmw[ind]
    ind = np.where(fmw > 0.8)[0]
    if ind.size > 0:
        axis_ratio_mh[ind] = (
            2.8-4.*axis_ratio_w[ind]+5*(axis_ratio_w[ind]-0.56)*fmw[ind])

    # make sure minimum axis ratio is enforced for small particles
    ind = np.where(diam < 6.)[0]
    if ind.size > 0:
        axis_ratio_mh[ind] = np.maximum(axis_ratio_w[ind], axis_ratio_mh[ind])

    return axis_ratio_mh


def compute_axis_ratio_dh(diam):
    """
    Computes the axis ratio of a dry hailstone

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)

    Returns
    -------
    axis_ratio : array of floats
       axis ratio

    """
    axis_ratio = np.zeros(diam.size)+0.8
    axis_ratio[diam < 10.] = 1.-0.02*diam[diam < 10.]

    return axis_ratio


def compute_axis_ratio_snow(diam):
    """
    Computes the axis ratio of snowflakes

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)

    Returns
    -------
    axis_ratio : array of floats
       axis ratio

    """
    axis_ratio = np.zeros(diam.size)+0.8

    return axis_ratio


def compute_axis_ratio_melting_snow(ar_ds, ar_rd, fmw):
    """
    Computes the axis ratio of melting snowflakes

    Parameters
    ----------
    ar_ds, ar_rd : array of floats
        axis ratio of dry snow and rain drops
    fmw : array of floats
        fraction of melted water

    Returns
    -------
    axis_ratio : array of floats
       axis ratio

    """
    return ar_ds+fmw*(ar_rd-ar_ds)


def compute_axis_ratio_ic(diam, coeff, exponent):
    """
    Computes the axis ratio of ice crystals

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)

    Returns
    -------
    axis_ratio : array of floats
       axis ratio

    """
    return coeff*np.power(diam, exponent-1.)


def compute_equi_vol_diam(mass, density=0.0009999720):
    """
    Given mass and density computes the equivalent volume diameter

    Parameters
    ----------
    mass : array of floats
        hydrometeor mass (g)
    density : array of floats
        hydrometeor density. The default is that of water at 4°C

    Returns
    -------
    diam : array of floats
       hydrometeor diameter

    """
    return np.power(6.*mass/(density*np.pi), 1./3.)


def compute_ms_core_diameter(dens_ms, dens_core, dens_shell, d_ms):
    """
    Computes the inner core diameter

    Parameters
    ----------
    dens_ms, dens_core, dens_shell : array of floats
        density of total melting snowflake, inner core and outer shell (g/mm3)
    d_ms : array of floats
        melting snowflake outer diameter

    Returns
    -------
    d_int : array of floats
       internal core diameter

    """
    alpha = np.power((dens_ms - dens_shell)/(dens_core-dens_shell), 1./3.)

    return d_ms*alpha


def compute_volume(diam):
    """
    Computes the volume of an spherical particle

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)

    Returns
    -------
    vol : array of floats
       particle volume (mm3)

    """
    return 4./3.*np.pi*np.power(diam/2., 3.)


def compute_ice_mass(diam, dens_ice=0.000916):
    """
    Given its equivalent diameter, computes the mass of an ice particle

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    dens_ice : float
        Density of the particle (g/mm3)

    Returns
    -------
    mass : array of floats
       mass of each particle (g)

    """
    return dens_ice*compute_volume(diam)


def compute_snow_mass(diam, dens_ice=0.000916):
    """
    Given its equivalent diameter, computes the mass of an aggregate

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    dens_ice : float
        Ice density (g/mm3)

    Returns
    -------
    mass : array of floats
       mass of each particle (g)

    """
    mass = np.zeros(diam.size)

    ind = np.where(diam < 2.)[0]
    if ind.size > 0:
        mass[ind] = 0.003*np.power(0.1*diam[ind], 2.)
    ind = np.where((diam >= 2.) & (diam < 20.))[0]
    if ind.size > 0:
        mass[ind] = 0.0067*np.power(0.1*diam[ind], 2.5)
    ind = np.where(diam >= 20.)[0]
    if ind.size > 0:
        mass[ind] = 0.0047*np.power(0.1*diam[ind], 3.)

    # constraint density range
    vol = compute_volume(diam)
    dens = mass/vol  # g/mm3

    ind = np.where(dens < 0.000005)[0]
    if ind.size > 0:
        mass[ind] = 0.000005*vol[ind]
    ind = np.where(dens > dens_ice)[0]
    if ind.size > 0:
        mass[ind] = dens_ice*vol[ind]

    return mass


def compute_hail_water_max(mass_ice, mass_ws):
    """
    computes the maximum water content of the hailstone coat before the onset
    of shedding

    Parameters
    ----------
    mass_ice : array of floats
        mass of ice of the hailstone (g)
    mass_ws : array of floats
        mass of soaked water (g)

    Returns
    -------
    mass_wc_max : array of floats
       maximum water content of the hail coat (g)

    """
    return 0.268+0.1389*(mass_ice+mass_ws)


def compute_ice_crystal_density(diam, coeff, exponent, dens_ice=0.000916):
    """
    Compute the ice crystal density. Assumes a power law with diameter

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    coeff, exponent : float
        the coefficient and exponent of the power law
    dens_ice : float
        Ice density (g/mm3)

    Returns
    -------
    dens : array of floats
       ice crystal density

    """
    dens = coeff*np.power(diam, exponent)
    dens[dens > dens_ice] = dens_ice

    return dens


def compute_snow_density(diam, mass, alpha=0.5):
    """
    Compute the total snow density, that of its inner core and that of the
    outer layer

    Parameters
    ----------
    diam : array of floats
        Equivalent volume diameter (mm)
    mass : array of floats
        snow mass (g)
    alpha : float
        relative position of transition between inner core and outer shell

    Returns
    -------
    dens_t, dens_core, dens_shell : array of floats
       Total snowflake density, inner core density and outer shell density

    """
    vol_t = compute_volume(diam)
    dens_t = mass/vol_t

    diam_core = diam*alpha
    mass_core = compute_snow_mass(diam_core)
    dens_core = mass_core/compute_volume(diam_core)
    dens_shell = (
        (dens_t-dens_core*np.power(alpha, 3.))/(1-np.power(alpha, 3.)))

    return dens_t, dens_core, dens_shell


def compute_melting_snow_density(dens_snow, fmw, dens_water=0.0009999720):
    """
    Computes the density of a melting snowflake

    Parameters
    ----------
    dens_snow : array of floats
        dry snowflake density
    fmw : array of floats
        fraction of melted water
    dens_water : float
        water density at 4 deg C (g/mm3)

    Returns
    -------
    dens : array of floats
       snowflake density

    """
    return dens_snow*dens_water/(fmw*dens_snow+(1-fmw)*dens_water)


def compute_hailstone_density(diam, dens_ice=0.000916, dens_min=0.0006,
                              diam_max=35., variable_density=True):
    """
    Computation of hailstone density

    Parameters
    ----------
    diam : array of floats
        hailstone diameter (mm)
    dens_ice : float
        density of pure ice
    dens_min : float
        minimum density of the hailstone
    diam_max : float
        diameter at which the density becomes maximum
    variable_density : bool
        If True the hailstone density increases with size
        If False the hailstone density is always equal to that of pure ice

    Returns
    -------
    dens_hs : array of floats
        hailstone density

    """
    if not variable_density:
        return dens_ice+np.zeros(diam.size)

    dens = dens_min+(dens_ice-dens_min)*np.log(1+diam)/np.log(1.+diam_max)
    dens[diam > diam_max] = dens_ice

    return dens


def compute_hs_air_vol(diam, mass_hail, dens_hail, dens_ice=0.000916):
    """
    Computation of the air volume of a hailstone

    Parameters
    ----------
    diam : array of floats
        hailstone diameter (mm)
    mass_hail : array of floats
        hailstone mass (g)
    dens_hail : array of floats
        hailstone density (g/mm3)
    dens_ice : float
        pure ice density

    Returns
    -------
    vol_air : array of floats
        hailstone air volume (mm3)

    """
    vol_air = np.zeros(diam.size)
    ind = np.where(~np.isclose(dens_hail, dens_ice))[0]
    if ind.size < 0:
        # this avoids numerical artifacts
        vol_air[ind] = compute_volume(diam[ind])-mass_hail[ind]/dens_ice
    return vol_air


def compute_velocity_rain(diam, rho0=1.22, rho=1.22):
    """
    Computation of raindrops terminal velocity

    Parameters
    ----------
    diam : array of floats
        raindrops diameter (mm)
    rho0 : float
        Default air density at sea level (kg/m3)
    rho : float
        Air density at the location of the rain drop (kg/m3)

    Returns
    -------
    vel : array of floats
        raindrop terminal velocity

    """
    return 3.78*np.power(diam, 0.67)*np.power(rho0/rho, 0.375+0.025*diam)


def compute_velocity_ic(length, mass, area_ratio, rho=1.22, visc=1.81e-5,
                        delta_0=8., c_0=0.35):
    """
    Computation of ice crystals terminal velocity

    Parameters
    ----------
    length : array of floats
        maximum dimension (length in columns, diameter of the base in plates)
        (mm)
    mass : array of floats
        ice crystal mass (g)
    area_ratio : array of floats
        ice crystal area ratio
    rho : float
        Air density at the location of the rain drop (kg/m3)
    visc : float
        dynamic viscosity of air (Kg/(m s))
    delta_0, c_0 : float
        Coefficients used in the computation of the Reynolds number

    Returns
    -------
    vel : array of floats
        raindrop terminal velocity

    """
    m_kg = mass/1000.

    # Best number
    x = (rho*8.*m_kg*g)/(visc*visc*np.pi*np.power(area_ratio, 0.5))

    # Reynolds number
    re = delta_0*delta_0/4.*np.power(
        np.sqrt(1.+4.*np.sqrt(x)/(delta_0*delta_0*np.sqrt(c_0)))-1., 2.)

    return visc*re/(rho*length)


def compute_velocity_snow(diam):
    """
    Computation of dry snowflakes terminal velocity

    Parameters
    ----------
    diam : array of floats
        snowflake diameter (mm)

    Returns
    -------
    vel : array of floats (m/s)
        snowflake terminal velocity

    """
    vel = np.zeros(diam.size)+1.3
    vel[diam < 10.] = 0.3+0.5*(np.log10(diam[diam < 10.])+1.)

    return vel


def compute_velocity_melting_snow(diam, mass, fmw, dens_air0=1.22,
                                  dens_air=1.22, dens_water=0.0009999720):
    """
    Computation of melting snowflake terminal velocity

    Parameters
    ----------
    diam : array of floats
        snowflake diameter (mm)
    mass : array of floats
        snowflake mass (g)
    fmw : array of floats
        fraction of water mass over the total snowflake mass
    dens_air0 : float
        air density at sea level (standard atmosphere) (kg/m3)
    dens_air : float
        air density at snowflake location (kg/m3)
    dens_water : float
        water density at 4 deg C (g/mm3)

    Returns
    -------
    vel : array of floats (m/s)
        melting snowflake terminal velocity

    """
    vel_ds = compute_velocity_snow(diam)
    vel_rd = compute_velocity_rain(
        np.power(6*mass/(dens_water*np.pi), 1./3.), rho0=dens_air0,
        rho=dens_air)
    c_coeff = 0.5*(vel_rd/vel_ds-1.)
    g_coeff = vel_rd/vel_ds-c_coeff*fmw-c_coeff*fmw*fmw

    return vel_rd/g_coeff


def compute_velocity_melting_hail(d_ext, mass, fmw, alpha, dens_air0=1.22,
                                  dens_air=1.22, visc=1.81e-5):
    """
    Computation of melting hailstone terminal velocity

    Parameters
    ----------
    d_ext : array of floats
        hailstone external diameter (mm)
    mass : array of floats
        hailstone mass (g)
    fmw : array of floats
        fraction of water mass over the total hailstone mass
    alpha : array of floats
        Ratio between the mass of soaked water and the mass of pure ice in
        the hailstone core. The value is NaN if the hailstone is not yet
        soaked
    dens_air0 : float
        air density at sea level (standard atmosphere) (kg/m3)
    dens_air : float
        air density at hailstone location (kg/m3)
    dens_water : float
        water density at 4 deg C (g/mm3)
    visc : float
        dynamic viscosity of air (Kg/(m s))

    Returns
    -------
    vel : array of floats (m/s)
        hailstone terminal velocity
    nre : array o floats
        Reynolds number

    """
    d_ext_m = d_ext/1000.

    nre = compute_reynolds_number(mass, dens_air=dens_air, visc=visc)

    mass_ice = mass*(1.-fmw)  # ice mass in grams
    mass_water = mass - mass_ice  # total water mass in grams
    # mass of soaked water (g):
    mass_w_soaked = np.zeros(d_ext.size)
    mass_w_soaked[np.isnan(alpha)] = mass_water[np.isnan(alpha)]
    mass_w_soaked[~np.isnan(alpha)] = (
        mass_ice[~np.isnan(alpha)]*alpha[~np.isnan(alpha)])

    mass_w_coat = mass_water - mass_w_soaked
    mass_ic_kg = (mass_ice+mass_w_soaked)/1000.

    # kynematic viscosity of air
    visc_kin = visc/dens_air

    # fall speed of dry hail stones
    vel_ds = visc_kin*nre/d_ext_m

    # fall speed of melting hailstones
    veq = deepcopy(vel_ds)
    ind = np.where(nre < 5000.)[0]
    if ind.size > 0:
        # equivalent to that of a drop of the same size
        # using Brandes et al. (2002) relationship
        veq[ind] = (
            np.power(dens_air0/dens_air, 0.5) *
            (-0.1021+4.932e3*d_ext_m[ind]-9.551e5*np.power(d_ext_m[ind], 2.)
             + 7.934e7*np.power(d_ext_m[ind], 3.)
             - 2.362e9*np.power(d_ext_m[ind], 4.)))

    ind = np.where((nre >= 5000.) & (nre < 25000.))[0]
    if ind.size > 0:
        veq[ind] = (
            (4.8e3+4.832e6*mass_ic_kg[ind])
            * visc_kin/d_ext_m[ind]*np.power(dens_air0/dens_air, 0.5))

    # critical water mass for onset of shedding [g]
    mass_water_max = compute_hail_water_max(mass_ice, mass_w_soaked)

    # this condition should never happen
    ind = np.where(mass_water > mass_water_max)[0]
    if ind.size > 0:
        mass_water[mass_water > mass_water_max] = (
            mass_water_max[mass_water > mass_water_max])

    # vel = vel_ds+mass_w_coat*(veq-vel_ds)/mass_water_max

    vel = vel_ds+fmw*(veq-vel_ds)

    return vel, nre


def compute_reynolds_number(mass, dens_air=1.22, visc=1.81e-5):
    """
    Computation of Reynolds number for each hailstone

    Parameters
    ----------
    mass : array of floats
        hailstone mass (g)
    dens_air : float
        air density at hailstone location (kg/m3)
    visc : float
        dynamic viscosity of air (kg /(m s)

    Returns
    -------
    nre : array of floats
        Reynolds number for each hailstone

    """
    mass_kg = mass/1000.

    # Best number
    nbe = 8.*mass_kg*g*dens_air/(np.pi*visc*visc)

    # Reynolds number
    nre = np.empty(mass.size)
    nre[nbe < 3.46e8] = 0.448*np.power(nbe[nbe < 3.46e8], 0.5536)
    nre[nbe >= 3.46e8] = np.power(nbe[nbe >= 3.46e8]/0.6, 0.5)

    return nre


def compute_reynolds_number_snow(diam, fmw, axis_ratio, vel, dens_air=1.22,
                                 visc=1.81e-5):
    """
    Computation of Reynolds number for each hailstone

    Parameters
    ----------
    diam : array of floats
        diameter of melting snowflake
    fmw : array of floats
        fraction of melted water
    axis_ratio : array of floats
        melting snowflake axis ratio
    vel : array of floats
        melting snowflake terminal velocity (m/s)
    dens_air : float
        air density at hailstone location (kg/m3)
    visc : float
        dynamic viscosity of air (Kg/(m s))

    Returns
    -------
    nre : array of floats
        Reynolds number for each snowflake
    cms : array of floats
        capacitance of the melting snowflake

    """
    d_m = diam/1000.
    p_ms = np.pi*d_m  # snowflake perimeter

    lamb = 0.8+0.2*fmw
    epsilon = 1e-6+np.sqrt(1.-np.power(axis_ratio, 2.))
    cs0 = d_m*epsilon/(2.*np.arcsin(epsilon))
    cms = lamb*cs0
    omega_ms = (
        np.pi*np.power(d_m/2., 2.)
        * (2.+np.pi*axis_ratio/epsilon
           * np.log((1.+epsilon)/(1.-epsilon))))  # surface area
    l_ms = omega_ms/p_ms  # length parameter

    return l_ms*vel*dens_air/visc, cms


def compute_snow_melting(temp, d_ds, d_ms, mass, axis_ratio, fmw, rel_h=100,
                         dens_air0=1.22, dens_air=1.22, visc=1.81e-5, temp0=0,
                         p_air0=1, p_air=1, t0k=273.15, w_water=18e-3,
                         delta_h=10.):
    """
    Computes the increase in water mass due to melting of snowflakes

    Parameters
    ----------
    temp : float
        ambient temperature (deg Celsius)
    d_ds : array of floats
        dry snowflake diameter
    d_ms : array of floats
        melting snowflake diameter
    mass : array of floats
        hailstone mass (g)
    axis_ratio : array of floats
        axis ratio
    fmw : array of floats
        fraction of water mass over the total hailstone mass
    rel_h : float
        relative humidity of air
    dens_air0 : float
        air density at sea level (standard atmosphere) (kg/m3)
    dens_air : float
        air density at hailstone location (kg/m3)
    visc : float
        dynamic viscosity of air (Kg/(m s))
    temp0 : float
        reference temperature (deg Celsius)
    p_air0 : float
        air pressure at the surface (hPa)
    p_air : float
        air pressure at the location of the hailstone (hPa)
    t0k : float
        0 deg Celsius value in Kelvin
    w_water : float
        molecular weight of water (Kg/mol)
    delta_h : float
        decrease in altitude of the particle (m)

    Returns
    -------
    delta_mw : array of floats
       Increase in water mass due to hailstone melting (g)
    vel : array of floats
        snowflakes velocity (m/s)

    """
    if temp < 0.:
        return 0.

    tk = temp+t0k
    tdiff = temp-temp0

    # thermal conductivity of air [J/(m s K)]
    ka = (2.381+0.0071*tdiff)*1e-2

    # diffusivity of water vapor in air [m2/s]
    dv = 2.11e-5*np.power(tk/t0k, 1.94)*(p_air0/p_air)

    # latent enthalpy of vaporization [J/Kg]
    lv = 2.499e6*np.power(t0k/tk, 0.167+3.67e-4*tk)

    # latent enthalpy of melting
    lm = 3.335e5*(1.+0.006*tdiff-3.14e-5*np.power(tdiff, 2.))

    nsc = visc/dens_air/dv  # schmidt number

    # saturation vapor pressure at 0 deg C [Pa]
    esat0 = 6.1094e2*np.exp(17.625*temp0/(243.04+temp0))

    # vapor density at temp t0 [Kg/m3]
    rhov0 = rel_h/100.*w_water/gas_constant*esat0/t0k

    # saturation vapor pressure at ambient temperature [Pa=Kg/(m s2)]
    esat = 6.1094e2*np.exp(17.625*temp/(243.04+temp))

    # ambient vapor density at temp far from the hail stone [Kg/m3]
    rhov_far = rel_h/100.*w_water/gas_constant*esat/tk

    vel = compute_velocity_melting_snow(
        d_ds, mass, fmw, dens_air0=dens_air0, dens_air=dens_air)

    nre, cms = compute_reynolds_number_snow(
        d_ms, fmw, axis_ratio, vel, dens_air=dens_air, visc=visc)

    x = np.power(nsc, 1./3.)*np.power(nre, 1./2.)
    f_coeff = np.zeros(x.size)  # ventilation coefficient of a snowflake
    ind = np.where(x <= 1.)[0]
    if ind.size > 0:
        f_coeff[ind] = 1.+0.14*x[ind]*x[ind]
    ind = np.where(x > 1.)[0]
    if ind.size > 0:
        f_coeff[ind] = 0.86+0.28*x[ind]

    delta_q = -4.*np.pi*f_coeff*cms/vel*(ka*tdiff+lv*dv*(rhov_far-rhov0))

    # gained water mass [g] in a delta_h m fall
    return -delta_q/lm*1000.*delta_h, vel


def compute_hail_melting(temp, diam_ext, diam_int, mass, fmw, alpha,
                         rel_h=100, dens_air0=1.22, dens_air=1.22,
                         visc=1.81e-5, temp0=0, p_air0=1, p_air=1, t0k=273.15,
                         w_water=18e-3, delta_h=10.):
    """
    Computes the increase in water mass due to melting of the hailstone

    Parameters
    ----------
    temp : float
        ambient temperature (deg Celsius)
    diam_ext : array of floats
        Total hailstone diameter (mm)
    diam_int : array of floats
        diameter of the internal ice core of the hailstone (mm)
    mass : array of floats
        hailstone mass (g)
    fmw : array of floats
        fraction of water mass over the total hailstone mass
    alpha : array of floats or NaNs
        Ratio between the mass of soaked water and the mass of pure ice in
        the hailstone core. The value is NaN if the hailstone is not yet
        soaked
    rel_h : float
        relative humidity of air
    dens_air0 : float
        air density at sea level (standard atmosphere) (kg/m3)
    dens_air : float
        air density at hailstone location (kg/m3)
    visc : float
        dynamic viscosity of air (Kg/(m s))
    temp0 : float
        reference temperature (deg Celsius)
    p_air0 : float
        air pressure at the surface (hPa)
    p_air : float
        air pressure at the location of the hailstone (hPa)
    t0k : float
        0 deg Celsius value in Kelvin
    w_water : float
        molecular weight of water (Kg/mol)
    delta_h : float
        decrease in altitude of the particle (m)

    Returns
    -------
    delta_mw : array of floats
       Increase in water mass due to hailstone melting (g)
    vel : array of floats
        hailstones velocity (m/s)

    """
    d_ext_m = diam_ext/1000.
    d_int_m = diam_int/1000.
    tk = temp+t0k
    tdiff = temp-temp0

    # thermal conductivity of air [J/(m s K)]
    ka = (2.381+0.0071*tdiff)*1e-2

    # thermal conductivity of water [J/(m s K)]
    kw = 0.568*np.exp(
        0.003473*tdiff-3.823e-5*np.power(tdiff, 2.)
        + 1.087e-6*np.power(tdiff, 3.))

    # diffusivity of water vapor in air [m2/s]
    dv = 2.11e-5*np.power(tk/t0k, 1.94)*(p_air0/p_air)

    # latent enthalpy of vaporization [J/Kg]
    lv = 2.499e6*np.power(t0k/tk, 0.167+3.67e-4*tk)

    # latent enthalpy of melting
    lm = 3.335e5*(1.+0.006*tdiff-3.14e-5*np.power(tdiff, 2.))

    # thermal diffusivity of air
    kka = 9.1018e-11*np.power(tk, 2.)+8.8197e-8*tk-1.0654e-5

    npr = visc/dens_air/kka  # prandtl number
    nsc = visc/dens_air/dv  # schmidt number

    # saturation vapor pressure at 0 deg C [Pa]
    esat0 = 6.1094e2*np.exp(17.625*temp0/(243.04+temp0))

    # vapor density at temp t0 [Kg/m3]
    rhov0 = rel_h/100.*w_water/gas_constant*esat0/t0k

    # saturation vapor pressure at ambient temperature [Pa=Kg/(m s2)]
    esat = 6.1094e2*np.exp(17.625*temp/(243.04+temp))

    # ambient vapor density at temp far from the hail stone [Kg/m3]
    rhov_far = rel_h/100.*w_water/gas_constant*esat/tk

    t_step = 0.001
    if temp < t_step:
        ta_vec = np.array([0.])
    else:
        ta_vec = np.arange(0., temp+t_step, t_step)

    # saturation water density at temperatures between 0 and temp
    esata = 6.1094e2*np.exp(17.625*ta_vec/(243.04+ta_vec))

    # ambient vapor density at the temp of the surface of the hail stone
    # [Kg/m3]
    rhov_close_vec = rel_h/100.*w_water/gas_constant*esata/(ta_vec+t0k)

    vel, nre = compute_velocity_melting_hail(
        diam_ext, mass, fmw, alpha, dens_air0=dens_air0, dens_air=dens_air,
        visc=visc)

    # thermal ventilation coefficient
    fh = 0.78+0.308*np.power(npr, 1./3.)*np.power(nre, 1./2.)

    # vapor ventilation coefficient
    fv = 0.78+0.308*np.power(nsc, 1./3.)*np.power(nre, 1./2.)

    delta_q = np.zeros(diam_ext.size)

    # nre < 250
    ind = np.where(nre < 250.)[0]
    if ind.size > 0:
        delta_q[ind] = (
            -(4.*np.pi*d_ext_m[ind]/vel[ind])
            * (ka*tdiff*fh[ind]+lv*dv*(rhov_far-rhov0)*fv[ind]))

    # nre between 250 and 3000
    ind = np.where((nre >= 250.) & (nre <= 3000.))[0]
    if ind.size > 0:
        delta_q[ind] = (
            -(2.*np.pi*d_ext_m[ind]/vel[ind])
            * (ka*tdiff*fh[ind]+lv*dv*(rhov_far-rhov0)*fv[ind]))

    # nre between 3000 and 6000
    ind = np.where((nre > 3000.) & (nre < 6000.))[0]
    if ind.size > 0:
        # no water layer yet
        ind2 = np.where(d_ext_m[ind] == d_int_m[ind])[0]
        if ind2.size > 0:
            ind2 = ind[ind2]
            delta_q[ind2] = (
                -(2.*np.pi*d_ext_m[ind2]/vel[ind2])
                * (ka*tdiff*fh[ind2]+lv*dv*(rhov_far-rhov0)*fv[ind2]))

        # water layer created
        ind2 = np.where(d_ext_m[ind] != d_int_m[ind])[0]
        if ind2.size > 0:
            ind2 = ind[ind2]

#            # tdiff_a = (T0-Ta) where T0 temperature at the surface of ice
#            # core and Ta temperature at the surface of the hailstone
#            # this is an aproximation that considers rhova = rhov0
#            tdiff_a = -(
#                ka*fh[ind2]*temp+lv*dv*fv[ind2]*(rhov_far-rhov0)
#                / (ka*fh[ind2]+kw*d_int_m[ind2]/(d_ext_m[ind2]-d_int_m[ind2])))
#
#            delta_q[ind2] = (
#                2.*np.pi*d_int_m[ind2]*d_ext_m[ind2]*kw*tdiff_a
#                / (vel[ind2]*(d_ext_m[ind2]-d_int_m[ind2])))

            # we try to find a solution to eq. B2C in Ryzhkov paper by finding
            # numerically the temperature that minimizes the function
            ta = np.zeros(ind2.size)
            rhov_close = np.zeros(ind2.size)
            for ind_ta, ii in enumerate(ind2):
                func = np.abs(
                    d_int_m[ii]*kw*(temp0-ta_vec)/(d_ext_m[ii]-d_int_m[ii])
                    + ka*(temp-ta_vec)*fh[ii]
                    + lv*dv*(rhov_far-rhov_close_vec)*fv[ii])
                ind3 = np.where(func == np.min(func))[0]
                ta[ind_ta] = ta_vec[ind3[0]]
                rhov_close[ind_ta] = rhov_close_vec[ind3[0]]

            delta_q[ind2] = (
                -(2.*np.pi*d_ext_m[ind2]/vel[ind2])
                * (ka*(temp-ta)*fh[ind2]
                   + lv*dv*(rhov_far-rhov_close)*fv[ind2]))

    # nre between 6000 and 20000
    ind = np.where((nre >= 6000.) & (nre < 20000.))[0]
    if ind.size > 0:
        delta_q[ind] = (
            -(0.76*np.pi*d_int_m[ind]*np.power(nre[ind], 1./2.)/vel[ind])
            * (np.power(npr, 1./3.)*ka*tdiff
               + np.power(nsc, 1./3.)*lv*dv*(rhov_far-rhov0)))

    # nre > 20000.
    ind = np.where(nre > 20000.)[0]
    if ind.size > 0:
        delta_q[ind] = (
            -(0.57+9.0e-6*nre[ind])
            * np.pi*d_int_m[ind]*np.power(nre[ind], 1./2.)/vel[ind]
            * (np.power(npr, 1./3.)*ka*tdiff
               + np.power(nsc, 1./3.)*lv*dv*(rhov_far-rhov0)))

    # forces negative change of enthalpy
    delta_q[delta_q > 0.] = 0.

    # gained water mass [g] in a delta_h m fall
    delta_mw = -delta_q/lm*1000.*delta_h
    delta_mw[(d_int_m == 0.) | (fmw == 1.)] = 0.

    return delta_mw, vel


def snowflake_change(delta_mw, mass_snow, d_ds, dens_ds, dens_core_init,
                     dens_shell_init, fmw, ar_ds, ar_rd):
    """
    Computation of the parameters of the snowflake resultant from melting

    Parameters
    ----------
    delta_mw : array of floats
        Increase of melting water mass due to melting (g)
    mass_snow : array of floats
        snowflake mass (g)
    d_ds : array of floats
        diameter of dry snowflake (mm)
    dens_ds, dens_core_init, dens_shell_init : array of floats
        total snowflake density when still not melted, inner core and outer
        shell (g/mm3)
    fmw : array of floats
        Fraction mass water
    ar_ds, ar_rd : array of floats
        axis ratio of a dry snow flake and axis ratio of a raindrop from
        equivalent size

    Returns
    -------
    d_ms, d_core : array of floats
        melting snowflake diameter and inner core diameter (mm)
    ar_ms : array of floats
        axis ratio of the melting snowflake
    dens_ms, dens_core, dens_shell : array of floats or None
        total density of the melting snowflake, inner core and outer shell
    fmw_out : array of floats
        New fraction of mass water

    """
    fmw_out = deepcopy(fmw)
    fmw_out += delta_mw/mass_snow
    fmw_out[fmw_out > 1.] = 1.

    dens_ms = compute_melting_snow_density(dens_ds, fmw_out)
    dens_core = compute_melting_snow_density(dens_core_init, fmw_out)
    dens_shell = compute_melting_snow_density(dens_shell_init, fmw_out)

    d_ms = np.power(dens_ds/dens_ms, 1./3.)*d_ds

    ar_ms = compute_axis_ratio_melting_snow(ar_ds, ar_rd, fmw_out)

    d_core = compute_ms_core_diameter(dens_ms, dens_core, dens_shell, d_ms)

    return d_ms, d_core, ar_ms, dens_ms, dens_core, dens_shell, fmw_out


def hailstone_change(delta_mw, vol_soakable, alpha, fmw, mass_hs,
                     dens_water=0.0009999720, dens_ice=0.000916):
    """
    Computation the parameters of the hail stone resultant from melting

    Parameters
    ----------
    delta_mw : array of floats
        Increase of melting water mass due to melting (g)
    vol_soakable : array of floats
        soakable volume in the hailstone core (mm3)
    alpha : array of floats or NaNs
        Ratio between the mass of soaked water and the mass of pure ice in
        the hailstone core. The value is NaN if the hailstone is not yet
        soaked
    fmw : array of floats
        Fraction mass water
    mass_hs : array of floats
        total hailstone mass (g)
    dens_water : float
        water density (g/mm3)
    dens_ice : float
        density of solid ice (g/mm3)

    Returns
    -------
    d_ext_out : array of floats
        Total hailstone equivalent volume diameter of the new hailstone (mm)
    d_int_out : array of floats
        Icy core equivalent volume diameter (mm)
    alpha_out : array of floats or None
        Ratio between the mass of soaked water and the mass of pure ice in
        the hailstone if the hailstone is soaked with water. NaN otherwise
    fmw_out : array of floats
        New fraction of mass water
    mass_hs_out : array of floats
        Mass of the new hailstone (g)
    mass_ws_out : array of floats
        Mass of the shed water (g)

    """
    n_hs = delta_mw.size

    # new water mass [g]
    mass_water = mass_hs*fmw+delta_mw
    mass_water[mass_water >= mass_hs] = mass_hs[mass_water >= mass_hs]
    fmw_out = mass_water/mass_hs

    # new water volume [mm3]
    vol_water = mass_water/dens_water

    # new ice mass [g]
    mass_ice = mass_hs - mass_water

    # new ice volume [mm3]
    vol_ice = mass_ice/dens_ice

    d_ext_out = np.zeros(n_hs)
    d_int_out = np.zeros(n_hs)
    alpha_out = deepcopy(alpha)
    mass_hs_out = deepcopy(mass_hs)
    mass_ws_out = np.zeros(n_hs)

    mass_wc = np.zeros(n_hs)
    mass_w_soaked = np.zeros(n_hs)
    mass_wc_max = np.zeros(n_hs)
    vol_ice_core = np.zeros(n_hs)

    # hail stone completely melted
    d_int_out[fmw_out == 1.] = 0
    d_ext_out[fmw_out == 1.] = np.power(
        6./np.pi*mass_hs[fmw_out == 1.]/dens_water, 1./3.)

    ind_not_soaked = np.where((fmw_out < 1.) & (np.isnan(alpha_out)))[0]
    if ind_not_soaked.size > 0:
        # hailstone not yet soaked in water
        ind_soaking = np.where(
            vol_water[ind_not_soaked] >= vol_soakable[ind_not_soaked])[0]
        if ind_soaking.size > 0:
            # hailstone reached soaked status
            ind_soaking = ind_not_soaked[ind_soaking]

            mass_w_soaked[ind_soaking] = vol_soakable[ind_soaking]*dens_water
            alpha_out[ind_soaking] = (
                mass_w_soaked[ind_soaking]/mass_ice[ind_soaking])
            mass_wc[ind_soaking] = (
                mass_water[ind_soaking] - mass_w_soaked[ind_soaking])
            mass_wc_max[ind_soaking] = compute_hail_water_max(
                mass_ice[ind_soaking], mass_w_soaked[ind_soaking])

            ind_shedding = np.where(
                mass_wc[ind_soaking] > mass_wc_max[ind_soaking])[0]
            if ind_shedding.size > 0:
                # shedding loss of hail mass
                ind_shedding = ind_soaking[ind_shedding]
                mass_ws_out[ind_shedding] = (
                    mass_wc[ind_shedding]-mass_wc_max[ind_shedding])
                mass_hs_out[ind_shedding] = (
                    mass_hs[ind_shedding]-mass_ws_out[ind_shedding])
                mass_wc[ind_shedding] = mass_wc_max[ind_shedding]
                mass_water[ind_shedding] = (
                    mass_w_soaked[ind_shedding]+mass_wc_max[ind_shedding])
                fmw_out[ind_shedding] = (
                    mass_water[ind_shedding]/mass_hs_out[ind_shedding])

            vol_ice_core[ind_soaking] = (
                mass_w_soaked[ind_soaking]/dens_water+vol_ice[ind_soaking])
            d_int_out[ind_soaking] = np.power(
                6./np.pi*vol_ice_core[ind_soaking], 1./3.)
            d_ext_out[ind_soaking] = np.power(
                6./np.pi*(
                    vol_ice_core[ind_soaking]+mass_wc[ind_soaking]/dens_water),
                1./3.)

        ind_not_soaking = np.where(
            vol_water[ind_not_soaked] < vol_soakable[ind_not_soaked])[0]
        if ind_not_soaking.size > 0:
            # hailstone not yet soaked
            ind_not_soaking = ind_not_soaked[ind_not_soaking]

            vol_ice_core[ind_not_soaking] = (
                vol_soakable[ind_not_soaking]+vol_ice[ind_not_soaking])
            d_int_out[ind_not_soaking] = np.power(
                6./np.pi*vol_ice_core[ind_not_soaking], 1./3.)
            d_ext_out[ind_not_soaking] = d_int_out[ind_not_soaking]

    ind_soaked = np.where((fmw_out < 1) & (~np.isnan(alpha_out)))[0]
    if ind_soaked.size > 0:
        # hailstone soaked
        mass_w_soaked[ind_soaked] = mass_ice[ind_soaked]*alpha_out[ind_soaked]
        mass_wc[ind_soaked] = (
            mass_water[ind_soaked] - mass_w_soaked[ind_soaked])
        mass_wc_max[ind_soaked] = compute_hail_water_max(
                mass_ice[ind_soaked], mass_w_soaked[ind_soaked])

        ind_shedding = np.where(
            mass_wc[ind_soaked] > mass_wc_max[ind_soaked])[0]
        # shedding loss of hail mass
        if ind_shedding.size > 0:
            ind_shedding = ind_soaked[ind_shedding]
            mass_ws_out[ind_shedding] = (
                mass_wc[ind_shedding]-mass_wc_max[ind_shedding])
            mass_hs_out[ind_shedding] = (
                mass_hs[ind_shedding]-mass_ws_out[ind_shedding])
            mass_wc[ind_shedding] = mass_wc_max[ind_shedding]
            mass_water[ind_shedding] = (
                mass_w_soaked[ind_shedding]+mass_wc_max[ind_shedding])
            fmw_out[ind_shedding] = (
                mass_water[ind_shedding]/mass_hs_out[ind_shedding])

        vol_ice_core[ind_soaked] = (
            mass_w_soaked[ind_soaked]/dens_water+vol_ice[ind_soaked])
        d_int_out[ind_soaked] = np.power(
            6./np.pi*vol_ice_core[ind_soaked], 1./3.)
        d_ext_out[ind_soaked] = np.power(
            6./np.pi*(
                vol_ice_core[ind_soaked]+mass_wc[ind_soaked]/dens_water),
            1./3.)

    return (
        d_ext_out, d_int_out, alpha_out, fmw_out, mass_hs_out, mass_ws_out)
