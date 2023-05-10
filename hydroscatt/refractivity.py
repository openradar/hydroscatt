"""
refractivity
============

computation of refractivity of hydrometeors

.. autosummary::
    :toctree: generated/

    wavelength_to_band
    refractive_index_water
    refractive_index_ice
    refractive_index_ice_particle
    refractive_index_ic_func
    refractive_index_melting_snow_core
    refractive_index_melting_snow_shell
    refractive_index_melting_hail_core
    refractive_index_mix

"""

from warnings import warn

import numpy as np
from scipy.constants import speed_of_light

from part_descrip import compute_volume, compute_ice_crystal_density


def wavelength_to_band(wavelength):
    """
    Given a wavelength get the frequency band

    Parameters
    ----------
    wavelength : float
        the wavelength (mm)

    Returns
    -------
    band : str
       frequency band

    """
    f_ghz = speed_of_light/(wavelength*1e-3)/1e9  # frequency in GHz

    if 2. <= f_ghz < 4.:
        return 'S'
    if 4. <= f_ghz < 8.:
        return 'C'
    if 8. <= f_ghz <= 12.:
        return 'X'

    warn('Unknown frequency band')

    return None


def refractive_index_water(wavelength, temp, sal=0.):
    """
    Dielectric constant of saline water according to Ellison (2005)

    Parameters
    ----------
    wavelength : float
        wavelength (mm). Valid range between 1 and 1000 GHz
    temp : float
        temperature (deg Celsius). Valid range between -20 and 30 deg Celsius
    sal : float
        salinity (promille). Valid range between 0 and 40 promille

    Returns
    -------
    m_water : float
       refractive index of water

    """
    f_ghz = speed_of_light/(wavelength*1e-3)/1e9  # frequency in GHz

    a1 = 0.46606917e-02
    a2 = -0.26087876e-04
    a3 = -0.63926782e-05
    a4 = 0.63000075e+01
    a5 = 0.26242021e-02
    a6 = -0.42984155e-02
    a7 = 0.34414691e-04
    a8 = 0.17667420
    a9 = -0.20491560e-03
    a10 = 0.58366888e+03
    a11 = 0.12634992e+03
    a12 = 0.69227972e-01
    a13 = 0.38957681e-03
    a14 = 0.30742330e+03
    a15 = 0.12634992e+03
    a16 = 0.37245044e+01
    a17 = 0.92609781e-02
    a18 = -0.26093754e-01

    temp2 = temp*temp
    temp3 = temp2*temp
    temp4 = temp3*temp
    sal2 = sal*sal
    alfa0 = (6.9431+3.2841*sal-0.099486*sal2)/(84.850+69.024*sal+sal2)
    alfa1 = 49.843-0.2276*sal+0.00198*sal2
    rtq = 1.+alfa0*(temp-15.)/(alfa1+temp)
    r15 = sal*(37.5109+5.45216*sal+1.4409e-02*sal2)/(1004.75+182.283*sal+sal2)
    sigma35 = (
        2.903602+8.607e-02*temp+4.738817e-04*temp2-2.991e-06*temp3
        + 4.3041e-09*temp4)
    sigma = sigma35*rtq*r15
    es0 = 87.85306
    es = es0*np.exp(-0.00456992*temp-a1*sal-a2*sal2-a3*sal*temp)
    e1 = a4*np.exp(-a5*temp-a6*sal-a7*sal*temp)
    einf = a16+a17*temp+a18*sal
    tau1 = (a8+a9*sal)*np.exp(a10/(temp+a11))
    tau2 = (a12+a13*sal)*np.exp(a14/(temp+a15))
    tp = 2.*np.pi/1000.
    delta1 = es-e1
    delta2 = e1-einf
    eps = (
        delta1/(1.-1j*tp*f_ghz*tau1)
        + delta2/(1.-1j*tp*f_ghz*tau2)+einf
        + 17.9751j*sigma/f_ghz)

    m_r = np.sqrt((np.sqrt(eps.real*eps.real+eps.imag*eps.imag)+eps.real)/2.)
    m_i = np.sqrt((np.sqrt(eps.real*eps.real+eps.imag*eps.imag)-eps.real)/2.)

    return m_r+1j*m_i


def refractive_index_ice(wavelength, temp):
    """
    Computes the refractive index of pure ice at a given temperature and
    wavelength

    Parameters
    ----------
    wavelength : float
        wavelength (mm)
    temp : float
        temperature (deg Celsius)

    Returns
    -------
    m_ice : float
       refractive index of ice

    """
    f_ghz = speed_of_light/(wavelength*1e-3)/1e9
    temp_k = temp+273.15

    e_1 = 3.1184+9.1e-4*(temp_k+273.15)

    theta = 300./temp_k-1.
    alpha = (0.00504+0.0062*theta)*np.exp(-22.1*theta)

    beta_m = (
        0.0207/temp_k
        * np.exp(335./temp_k)/np.power(np.exp(335./temp_k)-1, 2.)
        + 1.16e-11*np.power(f_ghz, 2.))
    delta_beta = np.exp(-9.963+0.0372*(temp_k-273.16))

    beta = beta_m+delta_beta

    e_2 = alpha/f_ghz+beta*f_ghz

    m_r = np.sqrt((np.sqrt(e_1*e_1+e_2*e_2)+e_1)/2.)
    m_i = np.sqrt((np.sqrt(e_1*e_1+e_2*e_2)-e_1)/2.)

    return m_r+1j*m_i


def refractive_index_ice_particle(m_ice, dens_particle, dens_ice=0.000916):
    """
    Computes the refractive index of an ice particle

    Parameters
    ----------
    m_ice : float
        refractive index of pure ice
    dens_particle : array of floats
        ice particle density
    dens_ice : float
        pure ice density

    Returns
    -------
    m_particle : array of floats
       refractive index of the ice particle

    """
    ks = dens_particle/dens_ice*(m_ice*m_ice-1.)/(m_ice*m_ice+2.)
    return np.sqrt((1.+2.*ks)/(1.-ks))


def refractive_index_ic_func(length, dens_coeff, dens_exponent, m_ice):
    """
    Computes the refractive index of an ice particle. Used to create a function
    that relates particle maximum dimension and refractive index

    Parameters
    ----------
    length : float
        maximum dimension (mm)
    dens_coeff, dens_exponent : float
        coefficient and exponent of the power law that relates maximum
        dimension with particle density
    m_ice : float
        pure ice refractive index

    Returns
    -------
    m_particle : float
       refractive index of the ice particle

    """
    dens_ic = compute_ice_crystal_density(
        np.array([length]), dens_coeff, dens_exponent)
    m_particle = refractive_index_ice_particle(m_ice, dens_ic)[0]

    return m_particle


def refractive_index_melting_snow_core(m_ice, m_water, m_air, dens_ds_core,
                                       dens_ms_core, fmw, d_core,
                                       dens_water=0.0009999720,
                                       dens_ice=0.0009167, dens_air=1.22e-6):
    """
    Computes the refractive index of the inner core of melting snowflakes

    Parameters
    ----------
    m_ice, m_water, m_air : float
        refractive index of pure ice, water and air
    dens_ds_core, dens_ms_core : array of float
        density of the dry snowflake and the melting snowflake inner core
    fmw : array of floats
        fraction mass water in the snowflake
    d_core : array of floats
        diameter of the internal core (mm)
    dens_water : float
        water density at 4 deg C g/mm3
    dens_ice : float
        ice density (g/mm3)
    dens_air : float
        air density (g/mm3)

    Returns
    -------
    m_core : array of floats
       refractive index of the snowflake inner core

    """
    fvw = dens_ds_core*fmw/(dens_water*(1-fmw)+dens_ds_core*fmw)
    vol_core = compute_volume(d_core)
    vol_water = fvw*vol_core
    vol_snow = (1.-fvw)*vol_core
    vol_ice = vol_snow*dens_ds_core/dens_ice
    fv_i_w = vol_ice/(vol_water+vol_ice)  # fraction volume of ice in water
    dens_iw = fvw*dens_water+(1-fvw)*dens_ice  # density of ice and water
    fv_a_iw = (dens_ms_core-dens_iw)/(dens_air-dens_iw)

    m_iw = refractive_index_mix(m_ice, m_water, fv_i_w)
    return refractive_index_mix(m_air, m_iw, fv_a_iw)


def refractive_index_melting_snow_shell(m_ice, m_water, m_air, dens_ds_shell,
                                        dens_ms_shell, fmw, d_core, d_shell,
                                        dens_water=0.0009999720,
                                        dens_ice=0.0009167, dens_air=1.22e-6):
    """
    Computes the refractive index of the outer shell of melting snowflakes

    Parameters
    ----------
    m_ice, m_water, m_air : float
        refractive index of pure ice, water and air
    dens_ds_shell, dens_ms_shell : array of float
        density of the dry snowflake and the melting snowflake outer shell
    fmw : array of floats
        fraction mass water in the snowflake
    d_core, d_shell : array of floats
        diameter of the internal core and external shell (mm)
    dens_water : float
        water density at 4 deg C g/mm3
    dens_ice : float
        ice density (g/mm3)
    dens_air : float
        air density (g/mm3)

    Returns
    -------
    m_core : array of floats
       refractive index of snowflake outer shell

    """
    fvw = dens_ds_shell*fmw/(dens_water*(1-fmw)+dens_ds_shell*fmw)
    vol_shell = compute_volume(d_shell)-compute_volume(d_core)
    vol_water = fvw*vol_shell
    vol_snow = (1.-fvw)*vol_shell
    vol_ice = vol_snow*dens_ds_shell/dens_ice
    fv_i_w = vol_ice/(vol_water+vol_ice)  # fraction volume of ice in water
    dens_iw = fvw*dens_water+(1-fvw)*dens_ice  # density of ice and water
    fv_iw_a = (dens_ms_shell-dens_air)/(dens_iw-dens_air)

    m_iw = refractive_index_mix(m_ice, m_water, fv_i_w)
    return refractive_index_mix(m_iw, m_air, fv_iw_a)


def refractive_index_melting_hail_core(wavelength, mass_hail, alpha, fmw,
                                       diam, dens_water=0.0009999720):
    """
    Computes the refractive index of the core of melting hailstones

    Parameters
    ----------
    wavelength : float
        wavelength (mm)
    mass_hail : array of floats
        total hailstone mass (g)
    alpha : array of floats
        ratio between the mass of soaked water and the mass of pure ice in
        the hailstone core. The value is NaN if the hailstone is not yet
        soaked
    fmw : array of floats
        fraction mass water in the hailstone
    diam : array of floats
        diameter of the internal ice core (mm)
    dens_water : float
        water density at 4 deg C g/mm3

    Returns
    -------
    m_core : array of floats
       refractive index of hailstone core

    """
    mass_ice_in_core = mass_hail*(1.-fmw)
    vol_core = np.pi/6*np.power(diam, 3.)

    # water mass in core
    mass_water_in_core = np.zeros(mass_hail.size)

    # core not yet soaked in water
    ind = np.where(np.isnan(alpha))[0]
    if ind.size > 0:
        mass_water_in_core[ind] = fmw[ind]*mass_hail[ind]

    # core soaked in water
    ind = np.where(~np.isnan(alpha))[0]
    if ind.size > 0:
        mass_water_in_core[ind] = alpha[ind]*mass_ice_in_core[ind]

    vol_water_in_core = mass_water_in_core/dens_water
    fvw_in_core = vol_water_in_core/vol_core
    dens_ice_in_core = mass_ice_in_core/(vol_core-vol_water_in_core)

    m_pure_ice = refractive_index_ice(wavelength, -1)
    m_ice_in_core = refractive_index_ice_particle(m_pure_ice, dens_ice_in_core)
    m_water_in_core = refractive_index_water(wavelength, 0.)

    return refractive_index_mix(m_water_in_core, m_ice_in_core, fvw_in_core)


def refractive_index_mix(m_1, m_2, fv1):
    """
    Computes the refractive index of a mixture of materials

    Parameters
    ----------
    m_1, m_2 : float
        refractive index of each material
    fv1 : array of floats
        fraction volume of first element with respect to the total

    Returns
    -------
    m_mix : array of floats
       resultant refractive index

    """
    er_1 = m_1.real*m_1.real-m_1.imag*m_1.imag
    ei_1 = 2.*m_1.real*m_1.imag
    e_1 = er_1+1j*ei_1

    er_2 = m_2.real*m_2.real-m_2.imag*m_2.imag
    ei_2 = 2.*m_2.real*m_2.imag
    e_2 = er_2+1j*ei_2

    e_12 = (
        e_2
        * (1.+2.*fv1*(e_1-e_2)/(e_1+2.*e_2))/(1.-fv1*(e_1-e_2)/(e_1+2.*e_2)))

    return np.sqrt(e_12)
