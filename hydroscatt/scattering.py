"""
scattering
============

computation of scattering properties of hydrometeors

.. autosummary::
    :toctree: generated/

    compute_tmatrix_2l
    compute_angular_moments
    compute_angular_moments_analytical
    compute_scattering_canting_sp
    compute_scattering_canting_psd
    compute_scattering_sp_tm
    compute_scattering_psd_tm
    compute_scattering_mixture

"""

import os
from shutil import copy

import numpy as np
import pandas as pd

from pytmatrix import tmatrix_aux, radar, scatter

from f_2l_tmatrix.tm2l import tmatrix as ftm2l_scattering


def compute_tmatrix_2l(fname_freq, fname_model_scatt, tm2l_dir, fname_scatt):
    """
    Using the fortran program tmatrix_py.f compute the scattering
    properties of a two-layer particle

    Parameters
    ----------
    fname_freq : str
        Name of the file containing the wavelength
    fname_model_scatt : str
        Name of the file containing the particle properties
    tm2l_dir : str
        Directory where the tmatrix_py.f program is placed
    fname_scatt : str
        Name of the file where to store the output

    Returns
    -------
    ang_moments_dict : dict
       dictionary containing the angular moments

    """
    copy(fname_model_scatt, f'{tm2l_dir}/fort.20')
    copy(fname_freq, f'{tm2l_dir}/spongyice.inp')
    cwd = os.getcwd()
    os.chdir(tm2l_dir)
    ftm2l_scattering(2)
    copy(f'{tm2l_dir}/spongyice.out', fname_scatt)
    os.remove(f'{tm2l_dir}/pck_tmatrix2b.out')
    os.remove(f'{tm2l_dir}/spongyice.out')
    os.remove(f'{tm2l_dir}/spongyice.inp')
    os.remove(f'{tm2l_dir}/fort.20')
    os.remove(f'{tm2l_dir}/fort.8')
    os.chdir(cwd)

    return fname_scatt


def compute_angular_moments(canting_angle, ele=0., n_part=500000, seed=0):
    """
    Compute the angular moments

    Parameters
    ----------
    canting_angle : array of floats
        the canting angle (deg)
    ele : float
        the elevation angle
    n_part : float
        number of particles used to compute the moments
    seed : int
        seed used to activate the random generator

    Returns
    -------
    ang_moments_dict : dict
       dictionary containing the angular moments

    """
    rng = np.random.default_rng(seed=seed)

    thet0 = np.pi/2-ele*np.pi/180.

    n_cant = canting_angle.size
    a1 = np.zeros(n_cant)
    a2 = np.zeros(n_cant)
    a3 = np.zeros(n_cant)
    a4 = np.zeros(n_cant)
    a5 = np.zeros(n_cant)
    a7 = np.zeros(n_cant)

    alpha = rng.uniform(
        low=0.0, high=360.0, size=n_part)*np.pi/180.
    for ind, cant_angle in enumerate(canting_angle):
        beta = rng.normal(
            scale=cant_angle, size=n_part)*np.pi/180.

        a1_aux = np.power(
            np.cos(beta)*np.sin(thet0)
            - np.sin(beta)*np.cos(thet0)*np.cos(alpha), 2.)
        a2_aux = np.sin(beta)*np.sin(beta)*np.sin(alpha)*np.sin(alpha)
        a3_aux = a1_aux*a1_aux
        a4_aux = a2_aux*a2_aux
        a5_aux = a1_aux*a2_aux
        a7_aux = a1_aux-a2_aux

        a1[ind] = np.mean(a1_aux)
        a2[ind] = np.mean(a2_aux)
        a3[ind] = np.mean(a3_aux)
        a4[ind] = np.mean(a4_aux)
        a5[ind] = np.mean(a5_aux)
        a7[ind] = np.mean(a7_aux)

    ang_moments_dict = {
        'a1': a1,
        'a2': a2,
        'a3': a3,
        'a4': a4,
        'a5': a5,
        'a7': a7
    }
    return ang_moments_dict


def compute_angular_moments_analytical(canting_angle):
    """
    Compute the angular moments using an analytical form. This form is only
    valid is the mean canting angle is equal to 0 and at low elevation angles

    Parameters
    ----------
    canting_angle : array of floats
        the canting angle (deg)

    Returns
    -------
    ang_moments_dict : dict
       dictionary containing the angular moments

    Reference
    ---------
    Ryzhkov et al., Polarimetric Radar Observation Operator for a Cloud Model
    with Spectral Microphysics, 2011, JAMC Vol. 50

    """
    cant_angl_rad = canting_angle*np.pi/180.
    r = np.exp(-2.*cant_angl_rad*cant_angl_rad)

    a1 = 0.25*np.power(1.+r, 2.)
    a2 = 0.25*(1.-np.power(r, 2.))
    a3 = np.power(3./8.+0.5*r+1./8.*np.power(r, 4.), 2.)
    a4 = (
        (3./8.-0.5*r+1./8.*np.power(r, 4.))
        * (3./8.+0.5*r+1./8.*np.power(r, 4.)))
    a5 = (
        1./8.*(3./8.+0.5*r+1./8.*np.power(r, 4.))
        * (1.-np.power(r, 4.)))
    a7 = 0.5*r*(1+r)

    ang_moments_dict = {
        'a1': a1,
        'a2': a2,
        'a3': a3,
        'a4': a4,
        'a5': a5,
        'a7': a7
    }
    return ang_moments_dict


def compute_scattering_canting_sp(wavelength, fv180, fh180, fv0, fh0,
                                  ang_moments_dict, kw_sqr=0.93,
                                  var_list=['sca_xsect_h', 'sca_xsect_v',
                                            'ext_xsect_h', 'ext_xsect_v',
                                            'refl_h', 'refl_v', 'ldr_h',
                                            'ldr_v', 'zdr', 'rho_hv',
                                            'delta_hv', 'kdp', 'A_h', 'A_v',
                                            'Adp']):
    """
    Compute scattering quantities of single particles averaged by canting
    angle

    Parameters
    ----------
    wavelength : float
        wave length (mm)
    fv180, fh180, fv0, fh0 : array of floats
        elements of the scattering matrix for each particle diameter
    ang_moments_dict : dict
        dictionary containing the angular moments
    kw_sqr : float
        the squared reference water dielectric factor for computing radar
        reflectivity
    var_list : list of str
        list of variables to compute

    Returns
    -------
    df : Pandas DataFrame
       a Pandas DataFrame containing all computed quantities for each particle
       diameter

    """
    a1 = ang_moments_dict['a1']
    a2 = ang_moments_dict['a2']
    a3 = ang_moments_dict['a3']
    a4 = ang_moments_dict['a4']
    a5 = ang_moments_dict['a5']
    a7 = ang_moments_dict['a7']

    if ('sca_xsect_h' in var_list or 'sca_xsect_v' in var_list
            or 'rho_hv' in var_list or 'ldr_h' in var_list
            or 'refl_h' in var_list or 'refl_v' in var_list
            or 'zdr' in var_list or 'delta_hv' in var_list
            or 'rho_hv' in var_list):
        j0 = np.power(np.abs(fh180), 2.)
        j1 = np.power(np.abs(fh180-fv180), 2.)
        j2 = np.conjugate(fh180)*(fh180-fv180)

    if ('ext_xsect_h' in var_list or 'ext_xsect_v' in var_list
            or 'A_h' in var_list or 'A_v' in var_list
            or 'Adp' in var_list or 'kdp' in var_list
            or 'zdr' in var_list):
        j3 = fh0-fv0
    if ('delta_hv' in var_list or 'rho_hv' in var_list):
        delta_co = j0+j1*a5-j2*a1-np.conjugate(j2)*a2

    if ('sca_xsect_h' in var_list or 'rho_hv' in var_list
            or 'ldr_h' in var_list or 'refl_h' in var_list
            or 'zdr' in var_list):
        sca_xsect_h = 4*np.pi*(j0-2.*j2.real*a2+j1*a4)
    if ('sca_xsect_v' in var_list or 'rho_hv' in var_list
            or 'ldr_v' in var_list or 'refl_v' in var_list
            or 'zdr' in var_list):
        sca_xsect_v = 4*np.pi*(j0-2.*j2.real*a1+j1*a3)
    if ('ext_xsect_h' in var_list or 'A_h' in var_list):
        ext_xsect_h = 2.*wavelength/1000.*(fh0.imag-j3.imag*a2)
    if ('ext_xsect_v' in var_list or 'A_v' in var_list):
        ext_xsect_v = 2.*wavelength/1000.*(fh0.imag-j3.imag*a3)
    if 'delta_hv' in var_list:
        delta_hv = -np.angle(delta_co, deg=True)
    if 'rho_hv' in var_list:
        rho_hv = 4.*np.pi*np.abs(delta_co)/np.sqrt(sca_xsect_h*sca_xsect_v)
    if 'kdp' in var_list:
        kdp = j3.real*a7*180./np.pi*wavelength
    if 'A_h' in var_list:
        # This is the two-way specific attenuation
        ah = 2*4.34e3*ext_xsect_h
    if 'A_v' in var_list:
        av = 2*4.34e3*ext_xsect_v
    if 'ldr_h' in var_list:
        ldr_h = 4*np.pi*j1*a5/sca_xsect_h
    if 'ldr_v' in var_list:
        ldr_v = 4*np.pi*j1*a5/sca_xsect_v
    if 'refl_h' in var_list:
        refl_h = (
            np.power(wavelength, 4.)*1e6/(np.power(np.pi, 5.)*kw_sqr)
            * sca_xsect_h)
    if 'refl_v' in var_list:
        refl_v = (
            np.power(wavelength, 4.)*1e6/(np.power(np.pi, 5.)*kw_sqr)
            * sca_xsect_v)
    if 'zdr' in var_list:
        zdr = sca_xsect_h/sca_xsect_v

    df = pd.DataFrame()
    if 'sca_xsect_h' in var_list:
        df['sca_xsect_h'] = 10.*np.log10(sca_xsect_h)
    if 'sca_xsect_v' in var_list:
        df['sca_xsect_v'] = 10.*np.log10(sca_xsect_v)
    if 'ext_xsect_h' in var_list:
        df['ext_xsect_h'] = 10.*np.log10(ext_xsect_h)
    if 'ext_xsect_v' in var_list:
        df['ext_xsect_v'] = 10.*np.log10(ext_xsect_v)
    if 'ldr_h' in var_list:
        df['ldr_h'] = 10.*np.log10(ldr_h)
    if 'ldr_v' in var_list:
        df['ldr_v'] = 10.*np.log10(ldr_v)
    if 'refl_h' in var_list:
        df['refl_h'] = 10.*np.log10(refl_h)
    if 'refl_v' in var_list:
        df['refl_v'] = 10.*np.log10(refl_v)
    if 'zdr' in var_list:
        df['zdr'] = 10.*np.log10(zdr)
    if 'rho_hv' in var_list:
        df['rho_hv'] = rho_hv
    if 'delta_hv' in var_list:
        df['delta_hv'] = delta_hv
    if 'kdp' in var_list:
        df['kdp'] = kdp
    if 'A_h' in var_list:
        df['A_h'] = ah
    if 'A_v' in var_list:
        df['A_v'] = av
    if 'Adp' in var_list:
        df['Adp'] = ah-av

    return df


def compute_scattering_canting_psd(wavelength, fv180, fh180, fv0, fh0,
                                   ang_moments_dict, delta_d, psd_vals,
                                   kw_sqr=0.93,
                                   var_list=['sca_xsect_h', 'sca_xsect_v',
                                             'ext_xsect_h', 'ext_xsect_v',
                                             'refl_h', 'refl_v', 'ldr_h',
                                             'ldr_v', 'zdr', 'rho_hv',
                                             'delta_hv', 'kdp', 'A_h', 'A_v',
                                             'Adp']):
    """
    Compute scattering quantities over a particle size distribution averaged
    by canting angle

    Parameters
    ----------
    wavelength : float
        wave length (mm)
    fv180, fh180, fv0, fh0 : array of floats
        elements of the scattering matrix for each particle diameter
    ang_moments_dict : dict
        dictionary containing the angular moments
    delta_d : array of floats
        particle diameter bin size
    psd_vals : array of floats
        Value of the particle size distribution for each particle diameter
    kw_sqr : float
        the squared reference water dielectric factor for computing radar
        reflectivity
    var_list : list of str
        list of variables to compute

    Returns
    -------
    scatt_dict : dict
       a dictionary containing all computed quantities

    """
    a1 = ang_moments_dict['a1']
    a2 = ang_moments_dict['a2']
    a3 = ang_moments_dict['a3']
    a4 = ang_moments_dict['a4']
    a5 = ang_moments_dict['a5']
    a7 = ang_moments_dict['a7']

    if ('sca_xsect_h' in var_list or 'sca_xsect_v' in var_list
            or 'rho_hv' in var_list or 'ldr_h' in var_list
            or 'ldr_v' in var_list or 'refl_h' in var_list
            or 'refl_v' in var_list or 'zdr' in var_list
            or 'delta_hv' in var_list or 'rho_hv' in var_list):
        j0 = np.power(np.abs(fh180), 2.)
        j1 = np.power(np.abs(fh180-fv180), 2.)
        j2 = np.conjugate(fh180)*(fh180-fv180)

    if ('ext_xsect_h' in var_list or 'ext_xsect_v' in var_list
            or 'A_h' in var_list or 'A_v' in var_list
            or 'Adp' in var_list or 'kdp' in var_list
            or 'zdr' in var_list):
        j3 = fh0-fv0
    if ('delta_hv' in var_list or 'rho_hv' in var_list):
        delta_co = j0+j1*a5-j2*a1-np.conjugate(j2)*a2

    if ('sca_xsect_h' in var_list or 'rho_hv' in var_list
            or 'ldr_h' in var_list or 'refl_h' in var_list
            or 'zdr' in var_list):
        sca_xsect_h = 4*np.pi*np.sum(
            (j0-2.*j2.real*a2+j1*a4)*psd_vals*delta_d)
    if ('sca_xsect_v' in var_list or 'rho_hv' in var_list
            or 'refl_v' in var_list or 'zdr' in var_list):
        sca_xsect_v = 4*np.pi*np.sum(
            (j0-2.*j2.real*a1+j1*a3)*psd_vals*delta_d)
    if ('ext_xsect_h' in var_list or 'A_h' in var_list):
        ext_xsect_h = 2.*wavelength/1000.*np.sum(
            (fh0.imag-j3.imag*a2)*psd_vals*delta_d)
    if ('ext_xsect_v' in var_list or 'A_v' in var_list):
        ext_xsect_v = 2.*wavelength/1000.*np.sum(
            (fh0.imag-j3.imag*a3)*psd_vals*delta_d)
    if 'delta_hv' in var_list:
        delta_hv = -np.angle(np.sum(delta_co*psd_vals*delta_d), deg=True)
    if 'rho_hv' in var_list:
        rho_hv = (
            4.*np.pi*np.abs(np.sum(delta_co*psd_vals*delta_d))
            / np.sqrt(sca_xsect_h*sca_xsect_v))
    if 'kdp' in var_list:
        kdp = np.sum(j3.real*a7*psd_vals*delta_d)*180./np.pi*wavelength
    if 'A_h' in var_list:
        # this is the two-way specific attenuation
        ah = 2*4.34e3*ext_xsect_h
    if 'A_v' in var_list:
        av = 2*4.34e3*ext_xsect_v
    if 'ldr_h' in var_list:
        ldr_h = 4*np.pi*np.sum(j1*a5*psd_vals*delta_d)/sca_xsect_h
    if 'ldr_v' in var_list:
        ldr_v = 4*np.pi*np.sum(j1*a5*psd_vals*delta_d)/sca_xsect_v
    if 'refl_h' in var_list:
        refl_h = (
            np.power(wavelength, 4.)*1e6/(np.power(np.pi, 5.)*kw_sqr)
            * sca_xsect_h)
    if 'refl_v' in var_list:
        refl_v = (
            np.power(wavelength, 4.)*1e6/(np.power(np.pi, 5.)*kw_sqr)
            * sca_xsect_v)
    if 'zdr' in var_list:
        zdr = sca_xsect_h/sca_xsect_v

    scatt_dict = {}
    if 'sca_xsect_h' in var_list:
        scatt_dict.update({'sca_xsect_h': 10.*np.log10(sca_xsect_h)})
    if 'sca_xsect_v' in var_list:
        scatt_dict.update({'sca_xsect_v': 10.*np.log10(sca_xsect_v)})
    if 'ext_xsect_h' in var_list:
        scatt_dict.update({'ext_xsect_h': 10.*np.log10(ext_xsect_h)})
    if 'ext_xsect_v' in var_list:
        scatt_dict.update({'ext_xsect_v': 10.*np.log10(ext_xsect_v)})
    if 'ldr_h' in var_list:
        scatt_dict.update({'ldr_h': 10.*np.log10(ldr_h)})
    if 'ldr_v' in var_list:
        scatt_dict.update({'ldr_v': 10.*np.log10(ldr_v)})
    if 'refl_h' in var_list:
        scatt_dict.update({'refl_h': 10.*np.log10(refl_h)})
    if 'refl_v' in var_list:
        scatt_dict.update({'refl_v': 10.*np.log10(refl_v)})
    if 'zdr' in var_list:
        scatt_dict.update({'zdr': 10.*np.log10(zdr)})
    if 'rho_hv' in var_list:
        scatt_dict.update({'rho_hv': rho_hv})
    if 'delta_hv' in var_list:
        scatt_dict.update({'delta_hv': delta_hv})
    if 'kdp' in var_list:
        scatt_dict.update({'kdp': kdp})
    if 'A_h' in var_list:
        scatt_dict.update({'A_h': ah})
    if 'A_v' in var_list:
        scatt_dict.update({'A_v': av})
    if 'Adp' in var_list:
        scatt_dict.update({'Adp': ah-av})

    return scatt_dict


def compute_scattering_sp_tm(scatterer,
                             geom_back=tmatrix_aux.geom_horiz_back,
                             geom_forw=tmatrix_aux.geom_horiz_forw,
                             var_list=['sca_xsect_h', 'sca_xsect_v', 'refl_h',
                                       'refl_v', 'ldr_h', 'ldr_v', 'zdr',
                                       'rho_hv', 'delta_hv', 'ext_xsect_h',
                                       'ext_xsect_v', 'kdp', 'A_h', 'A_v',
                                       'Adp']):
    """
    Compute single particle scattering quantities using pytmatrix

    Parameters
    ----------
    scatterer : pytmatrix scatterer object
        object containing the scatterers information
    geom_back, geom_forw : tuple
        forward and backward scattering geometry
    var_list : list of str
        list of variables to compute

    Returns
    -------
    scatt_dict : dict
       a dictionary containing all computed quantities

    """
    scatt_dict = {}
    scatterer.set_geometry(geom_back)
    if 'sca_xsect_h' in var_list:
        scatt_dict.update({
            'sca_xsect_h': 10.*np.log10(
                scatter.sca_xsect(scatterer, h_pol=True)/1e6)})
    if 'sca_xsect_v' in var_list:
        scatt_dict.update({
            'sca_xsect_v': 10.*np.log10(
                scatter.sca_xsect(scatterer, h_pol=False)/1e6)})
    if 'refl_h' in var_list:
        scatt_dict.update({'refl_h': 10.*np.log10(radar.refl(scatterer))})
    if 'refl_v' in var_list:
        scatt_dict.update({
            'refl_v': 10.*np.log10(radar.refl(scatterer, False))})
    if 'zdr' in var_list:
        scatt_dict.update({'zdr': 10.*np.log10(radar.Zdr(scatterer))})
    if 'ldr_h' in var_list:
        scatt_dict.update({
            'ldr_h': 10.*np.log10(scatter.ldr(scatterer, h_pol=True))})
    if 'ldr_v' in var_list:
        scatt_dict.update({
            'ldr_v': 10.*np.log10(scatter.ldr(scatterer, False))})
    if 'rho_hv' in var_list:
        scatt_dict.update({'rho_hv': radar.rho_hv(scatterer)})
    if 'delta_hv' in var_list:
        scatt_dict.update({'delta_hv': 180./np.pi*radar.delta_hv(scatterer)})

    scatterer.set_geometry(geom_forw)
    if 'ext_xsect_h' in var_list:
        scatt_dict.update({
            'ext_xsect_h': 10.*np.log10(
                scatter.ext_xsect(scatterer, h_pol=True)/1e6)})
    if 'ext_xsect_v' in var_list:
        scatt_dict.update({
            'ext_xsect_v': 10.*np.log10(
                scatter.ext_xsect(scatterer, h_pol=False)/1e6)})
    if 'kdp' in var_list:
        scatt_dict.update({'kdp': radar.Kdp(scatterer)})
    if 'A_h' in var_list or 'Adp' in var_list:
        Ah = 2*radar.Ai(scatterer)
        if 'A_h' in var_list:
            scatt_dict.update({'A_h': Ah})
    if 'A_v' in var_list or 'Adp' in var_list:
        Av = 2*radar.Ai(scatterer, False)
        if 'A_v' in var_list:
            scatt_dict.update({'A_v': Av})
    if 'Adp' in var_list:
        scatt_dict.update({'Adp': Ah-Av})

    return scatt_dict


def compute_scattering_psd_tm(scatterer,
                              geom_back=tmatrix_aux.geom_horiz_back,
                              geom_forw=tmatrix_aux.geom_horiz_forw,
                              var_list=['refl_h', 'refl_v', 'ldr_h', 'ldr_v',
                                        'zdr', 'rho_hv', 'delta_hv', 'kdp',
                                        'A_h', 'A_v', 'Adp']):
    """
    Compute PSD scattering quantities using pytmatrix

    Parameters
    ----------
    scatterer : pytmatrix scatterer object
        object containing the scatterers information
    geom_back, geom_forw : tuple
        forward and backward scattering geometry
    var_list : list of str
        list of variables to compute

    Returns
    -------
    scatt_dict : dict
       a dictionary containing all computed quantities

    """
    scatt_dict = {}
    scatterer.set_geometry(geom_back)
    if 'refl_h' in var_list:
        scatt_dict.update({'refl_h': 10.*np.log10(radar.refl(scatterer))})
    if 'refl_v' in var_list:
        scatt_dict.update({
            'refl_v': 10.*np.log10(radar.refl(scatterer, False))})
    if 'zdr' in var_list:
        scatt_dict.update({'zdr': 10.*np.log10(radar.Zdr(scatterer))})
    if 'ldr_h' in var_list:
        scatt_dict.update({'ldr_h': 10.*np.log10(radar.ldr(scatterer))})
    if 'ldr_v' in var_list:
        scatt_dict.update({
            'ldr_v': 10.*np.log10(radar.ldr(scatterer, False))})
    if 'rho_hv' in var_list:
        scatt_dict.update({'rho_hv': radar.rho_hv(scatterer)})
    if 'delta_hv' in var_list:
        scatt_dict.update({'delta_hv': 180./np.pi*radar.delta_hv(scatterer)})

    scatterer.set_geometry(geom_forw)
    if 'kdp' in var_list:
        scatt_dict.update({'kdp': radar.Kdp(scatterer)})
    if 'A_h' in var_list or 'Adp' in var_list:
        Ah = 2*radar.Ai(scatterer)
        if 'A_h' in var_list:
            scatt_dict.update({'A_h': Ah})
    if 'A_v' in var_list or 'Adp' in var_list:
        Av = 2*radar.Ai(scatterer, False)
        if 'A_v' in var_list:
            scatt_dict.update({'A_v': Av})
    if 'Adp' in var_list:
        scatt_dict.update({'Adp': Ah-Av})

    return scatt_dict


def compute_scattering_mixture(scatt_dict_1, scatt_dict_2,
                               var_list=['refl_h', 'refl_v', 'ldr_h',
                                         'zdr', 'rho_hv', 'delta_hv', 'kdp',
                                         'A_h', 'A_v', 'Adp']):
    """
    Given the scattering properties of two hydrometeor species computes the
    properties of the mixture

    Parameters
    ----------
    scatt_dict_1, scatt_dict_2 : dict
        dictionary containing the scattering parameters of each hydrometeor
        species
    var_list : list of str
        list of variables to compute

    Returns
    -------
    scatt_dict : dict
       a dictionary containing all computed quantities

    """
    scatt_dict = {}
    if ('refl_h' in var_list or 'zdr' in var_list or 'rho_hv' in var_list
            or 'ldr_h' in var_list):
        refl_h = 10.*np.log10(
            np.power(10., 0.1*scatt_dict_1['refl_h'])
            + np.power(10., 0.1*scatt_dict_2['refl_h']))
        if 'refl_h' in var_list:
            scatt_dict.update({'refl_h': refl_h})
    if ('refl_v' in var_list or 'zdr' in var_list or 'rho_hv' in var_list
            or 'ldr_v' in var_list):
        refl_v = 10.*np.log10(
            np.power(10., 0.1*scatt_dict_1['refl_v'])
            + np.power(10., 0.1*scatt_dict_2['refl_v']))
        if 'refl_v' in var_list:
            scatt_dict.update({'refl_v': refl_v})
    if 'zdr' in var_list:
        scatt_dict.update({'zdr': refl_h-refl_v})
    if 'A_h' in var_list:
        scatt_dict.update({'A_h': scatt_dict_1['A_h']+scatt_dict_2['A_h']})
    if 'A_v' in var_list:
        scatt_dict.update({'A_v': scatt_dict_1['A_v']+scatt_dict_2['A_v']})
    if 'Adp' in var_list:
        scatt_dict.update({'Adp': scatt_dict_1['Adp']+scatt_dict_2['Adp']})
    if 'kdp' in var_list:
        scatt_dict.update({'kdp': scatt_dict_1['kdp']+scatt_dict_2['kdp']})
    if 'delta_hv' in var_list:
        scatt_dict.update({
            'delta_hv': scatt_dict_1['delta_hv']+scatt_dict_2['delta_hv']})
    if 'rho_hv' in var_list:
        rho1_nom = (
            scatt_dict_1['rho_hv']
            * np.sqrt(
                np.power(10., 0.1*scatt_dict_1['refl_h'])
                * np.power(10., 0.1*scatt_dict_1['refl_v'])))
        rho2_nom = (
            scatt_dict_2['rho_hv']
            * np.sqrt(
                np.power(10., 0.1*scatt_dict_2['refl_h'])
                * np.power(10., 0.1*scatt_dict_2['refl_v'])))
        scatt_dict.update({
            'rho_hv': (
                (rho1_nom+rho2_nom)
                / np.sqrt(
                    np.power(10., 0.1*refl_h)*np.power(10., 0.1*refl_v)))})
    if 'ldr_h' in var_list:
        ldr1_nom = (
            np.power(10., 0.1*scatt_dict_1['ldr_h'])
            * np.power(10., 0.1*scatt_dict_1['refl_h']))
        ldr2_nom = (
            np.power(10., 0.1*scatt_dict_2['ldr_h'])
            * np.power(10., 0.1*scatt_dict_2['refl_h']))
        scatt_dict.update({'ldr_h': 10.*np.log10(ldr1_nom+ldr2_nom) - refl_h})
    if 'ldr_v' in var_list:
        ldr1_nom = (
            np.power(10., 0.1*scatt_dict_1['ldr_v'])
            * np.power(10., 0.1*scatt_dict_1['refl_v']))
        ldr2_nom = (
            np.power(10., 0.1*scatt_dict_2['ldr_v'])
            * np.power(10., 0.1*scatt_dict_2['refl_v']))
        scatt_dict.update({'ldr_v': 10.*np.log10(ldr1_nom+ldr2_nom) - refl_v})

    return scatt_dict
