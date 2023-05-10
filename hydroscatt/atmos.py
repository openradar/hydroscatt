"""
atmos
=====

computation of parameters related to the atmosphere

.. autosummary::
    :toctree: generated/

    compute_air_pressure
    compute_air_density
    compute_dynamic_viscosity

"""

import numpy as np

# g : gravitational acceleration 9.80665 m/s2
# gas_constant: 8.314462618 J/(mol K)
from scipy.constants import g, gas_constant


def compute_air_pressure(alt, lapse_rate=6.5, p_air0=1.075, m_air=0.0289644,
                         alt0=-610, temp0=19., tok=273.15):
    """
    Computes the air pressure at a particular altitude above sea level.
    It uses the barometric formula. The default values are those of the
    International Standard Atmosphere

    Parameters
    ----------
    alt : float
        altitude (masl)
    lapse_rate : float
        temperature lapse rate (deg C/km)
    p_air0 : float
        reference air pressure (hPa)
    m_air : float
        molar mass of Earth's air (Kg/mol)
    alt0 : float
        reference altitude (masl)
    temp0 : float
        reference temperature (deg C)
    tok : float
        Conversion from deg C to Kelvin

    Returns
    -------
    p_air : float
       air pressure (hPa)

    References
    ----------
    Barometric formula:
        https://en.wikipedia.org/wiki/Barometric_formula
    International Standard Atmosphere:
        https://en.wikipedia.org/wiki/International_Standard_Atmosphere

    """
    temp0_k = temp0+tok
    lapse_rate_m = lapse_rate/1e3
    return p_air0*np.power(
        (temp0_k+(alt-alt0)*lapse_rate_m)/temp0_k,
        -g*m_air/(gas_constant*lapse_rate_m))


def compute_air_density(alt, lapse_rate=6.5, dens_air0=1.2985, m_air=0.0289644,
                        alt0=-610, temp0=19., tok=273.15):
    """
    Computes the air density at a particular altitude above sea level.
    It uses the barometric formula. The default values are those of the
    International Standard Atmosphere

    Parameters
    ----------
    alt : float
        altitude (masl)
    lapse_rate : float
        temperature lapse rate (deg C/km)
    dens_air0 : float
        reference air density (Kg/m3)
    m_air : float
        molar mass of Earth's air (Kg/mol)
    alt0 : float
        reference altitude (masl)
    temp0 : float
        reference temperature (deg C)
    tok : float
        Conversion from deg C to Kelvin

    Returns
    -------
    p_air : float
       air pressure (hPa)

    References
    ----------
    Barometric formula:
        https://en.wikipedia.org/wiki/Barometric_formula
    International Standard Atmosphere:
        https://en.wikipedia.org/wiki/International_Standard_Atmosphere

    """
    temp0_k = temp0+tok
    lapse_rate_m = lapse_rate/1e3
    return dens_air0*np.power(
        temp0_k/(temp0_k+(alt-alt0)*lapse_rate_m),
        1.+g*m_air/(gas_constant*lapse_rate_m))


def compute_dynamic_viscosity(temp, tok=273.15):
    """
    Computes the dynamic viscosity

    Parameters
    ----------
    temp : float
        temperature (deg C)
    tok : float
        Conversion from deg C to Kelvin

    Returns
    -------
    visc : float
       dynamic viscosity of air (Kg/(m s))

    References
    ----------
    https://en.wikipedia.org/wiki/Viscosity#Dynamic_(shear)_viscosity

    """
    temp_k = temp+tok
    return 2.791e-7*np.power(temp_k, 0.7355)
