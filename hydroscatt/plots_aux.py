#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
plots_aux
================================================

auxiliary plots

.. autosummary::
    :toctree: generated/

    plot_air_density
    plot_air_pressure
    plot_axis_ratio_rain
    psd_hail_aloft
    hail_density
    plot_snow_mass
    plot_equi_vol_diam
    plot_snow_density

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np

from pytmatrix.psd import ExponentialPSD
from pytmatrix.tmatrix_aux import dsr_thurai_2007

from part_descrip import compute_hailstone_density, compute_axis_ratio_brandes
from part_descrip import compute_snow_mass, compute_snow_density
from part_descrip import compute_equi_vol_diam
from atmos import compute_air_density, compute_air_pressure
from atmos import compute_dynamic_viscosity
from graph import plot_multiple_var, plot_vertical_profile
from scattering_io import get_save_dir

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to hail modelling framework')

    # keyword arguments
    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='path where to store the plots')

    args = parser.parse_args()

    print(f'====== auxiliary plots started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== auxiliary plots finished: ")

    savedir = get_save_dir(
        args.path, None, None, None, create_dir=True,
        with_subdirs=False)

    plot_dynamic_viscosity(savedir)
#    plot_axis_ratio_rain(savedir)
#    plot_snow_mass(savedir)
#    plot_equi_vol_diam(savedir)
#    plot_snow_density(savedir, logy=True)


def plot_air_density(savedir, alt_min=0., alt_max=4000., alt_step=10.):
    """
    Plots the air density as a function of altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    alt_min, alt_max, alt_step : float
        minimum and maximum altitude and step to plot (m)

    Returns
    -------
    Nothing

    """
    alt_vec = np.arange(alt_min, alt_max+alt_step, alt_step)
    dens_air = compute_air_density(alt_vec)
    fname = f'{savedir}air_density.png'
    plot_vertical_profile(
        [dens_air], alt_vec, xlabel='air density (kg/m3)',
        ylabel='altitude (masl)',
        titl='standard atmosphere air density', fname=fname,
        invert_yaxis=False)


def plot_air_pressure(savedir, alt_min=0., alt_max=4000., alt_step=10.):
    """
    Plots the air pressure as a function of altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    alt_min, alt_max, alt_step : float
        minimum and maximum altitude and step to plot (m)

    Returns
    -------
    Nothing

    """
    alt_vec = np.arange(alt_min, alt_max+alt_step, alt_step)
    p_air = compute_air_pressure(alt_vec)
    fname = f'{savedir}air_pressure.png'
    plot_vertical_profile(
        [p_air], alt_vec, xlabel='air pressure (hPa)',
        ylabel='altitude (masl)',
        titl='standard atmosphere air pressure', fname=fname,
        invert_yaxis=False)


def plot_dynamic_viscosity(savedir, alt_min=0., alt_max=4000., alt_step=10.,
                           lapse_rate=6.5, alt_iso0=4000.,):
    """
    Plots the air density as a function of altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    alt_min, alt_max, alt_step : float
        minimum and maximum altitude and step to plot (m)
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)

    Returns
    -------
    Nothing

    """
    alt_vec = np.arange(alt_min, alt_max+alt_step, alt_step)
    temp_vec = (alt_iso0-alt_vec)*lapse_rate/1000.
    visc = compute_dynamic_viscosity(temp_vec)
    fname = f'{savedir}air_viscosity.png'
    plot_vertical_profile(
        [visc], alt_vec, xlabel='dynamic viscosity of air (Kg / (m s))',
        ylabel='altitude (masl)',
        titl='standard atmosphere dynamic viscosity of air', fname=fname,
        invert_yaxis=False)


def plot_axis_ratio_rain(savedir, d_min=0.1, d_max=10, d_step=0.1,
                         relations=['Brandes', 'Thurai']):
    """
    Plots the axis ratio of raindrops

    Parameters
    ----------
    savedir : str
        directory where to save the data
    d_min, d_max, d_step: float
        minimum and maximum diameter and step to plot (mm)
    relations : list of str
        list of relationships to plot

    Returns
    -------
    Nothing

    """
    diam = np.arange(d_min, d_max+d_step, d_step)
    ar_list = []
    for relation in relations:
        if relation == 'Brandes':
            ar_list.append(compute_axis_ratio_brandes(diam))
        elif relation == 'Thurai':
            axis_ratio = np.zeros(diam.size)
            for ind, diam_part in enumerate(diam):
                axis_ratio[ind] = dsr_thurai_2007(diam_part)
            ar_list.append(axis_ratio)

    fname = f'{savedir}ar_rain.png'
    plot_multiple_var(
        [diam], ar_list, xlabel='diameter (mm)',
        ylabel='aspect ratio', titl='aspect ratio of rain',
        fname=fname, labels=relations)


def psd_hail_aloft(savedir, d_min=0.1, d_step=0.1):
    """
    Plots PSD hail aloft

    Parameters
    ----------
    savedir : str
        directory where to save the data
    d_min, d_step: float
        minimum diameter and step to plot

    Returns
    -------
    Nothing

    """
    # PSD graupel
    n_g = 8000.
    lamb_g = 1.6
    d_max_g = 8
    psd_g = ExponentialPSD(N0=n_g, Lambda=lamb_g, D_max=d_max_g)

    # PSD hail
    # corresponds to no hail aloft, small hail, moderate hail and large hail
    label_list = ['no hail', 'small hail', 'moderate hail', 'large hail']
    d_max_list = np.array([0., 14., 24., 35.])
    lamb_h_list = np.array([0., 0.99, 0.42, 0.27])
    coeff_nh_list = np.array([0., 200., 400., 800.])

    diam = np.arange(d_min, d_max_list.max()+d_step, d_step)
    psd_h_vals_list = []
    for d_max_h, lamb_h, coeff_nh in zip(d_max_list, lamb_h_list,
                                         coeff_nh_list):
        n_h = coeff_nh*np.power(lamb_h, 4.11)
        psd_h = ExponentialPSD(N0=n_h, Lambda=lamb_h, D_max=d_max_h)
        psd_h_vals = psd_h(diam) + psd_g(diam)
        psd_h_vals[psd_h_vals == 0.] = np.nan
        psd_h_vals_list.append(psd_h_vals)

    fname = f'{savedir}psd_hail_aloft.png'
    plot_multiple_var(
        [diam], psd_h_vals_list, xlabel='diameter (mm)',
        ylabel='N(D)(m-3,mm-1)', titl='hail size distribution aloft',
        fname=fname, labels=label_list, logy=True)


def hail_density(savedir, d_min=0.1, d_max=80, d_step=0.1):
    """
    Plots hail density aloft

    Parameters
    ----------
    savedir : str
        directory where to save the data
    d_min, d_max, d_step: float
        minimum and maximum diameter and step to plot

    Returns
    -------
    Nothing

    """
    label_list = ['fixed dens', 'var dens']
    diam = np.arange(d_min, d_max+d_step, d_step)
    fixed_dens = compute_hailstone_density(
        diam, dens_ice=0.000916, dens_min=0.0006, diam_max=35.,
        variable_density=False)

    var_dens = compute_hailstone_density(
        diam, dens_ice=0.000916, dens_min=0.0006, diam_max=35.,
        variable_density=True)

    fname = f'{savedir}hail_density_aloft.png'
    plot_multiple_var(
        [diam], [fixed_dens, var_dens], xlabel='diameter (mm)',
        ylabel='density (g/cm3)', titl='hail density aloft',
        fname=fname, labels=label_list)


def plot_snow_mass(savedir, d_min=0.1, d_max=20, d_step=0.1):
    """
    Plots snow mass aloft

    Parameters
    ----------
    savedir : str
        directory where to save the data
    d_min, d_max, d_step: float
        minimum and maximum diameter and step to plot

    Returns
    -------
    Nothing

    """
    diam = np.arange(d_min, d_max+d_step, d_step)
    mass_snow = compute_snow_mass(diam)

    fname = f'{savedir}snow_mass_aloft.png'
    plot_multiple_var(
        [diam], [mass_snow], xlabel='diameter (mm)',
        ylabel='mass (g)', titl='snow mass aloft',
        fname=fname, labels=None)


def plot_equi_vol_diam(savedir, d_min=0.1, d_max=20, d_step=0.1):
    """
    Plots hail density aloft

    Parameters
    ----------
    savedir : str
        directory where to save the data
    d_min, d_max, d_step: float
        minimum and maximum diameter and step to plot

    Returns
    -------
    Nothing

    """
    diam = np.arange(d_min, d_max+d_step, d_step)
    mass_snow = compute_snow_mass(diam)
    d_rd = compute_equi_vol_diam(mass_snow)

    fname = f'{savedir}snow_equi_vol_diam.png'
    plot_multiple_var(
        [diam], [d_rd], xlabel='snowflake diameter (mm)',
        ylabel='raindrop diameter (mm)',
        titl='equivalent raindrop diameter of a snowflake',
        fname=fname, labels=None)


def plot_snow_density(savedir, d_min=0.1, d_max=20, d_step=0.1, alpha=0.5,
                      logy=False):
    """
    Plots snow density aloft

    Parameters
    ----------
    savedir : str
        directory where to save the data
    d_min, d_max, d_step: float
        minimum and maximum diameter and step to plot
    alpha : float
        ratio between the internal core diameter and the external shell
        diameter
    logy : bool
        if True the y axis will be expressed in log

    Returns
    -------
    Nothing

    """
    label_list = ['snowflake', 'core', 'shell']
    ylabel = 'density (g/mm3)'
    if logy:
        ylabel = 'log density (g/mm3)'

    diam = np.arange(d_min, d_max+d_step, d_step)
    mass_snow = compute_snow_mass(diam)
    dens_ds, dens_ds_core, dens_ds_shell = compute_snow_density(
        diam, mass_snow, alpha=alpha)

    fname = f'{savedir}snow_density_aloft.png'
    plot_multiple_var(
        [diam], [dens_ds, dens_ds_core, dens_ds_shell],
        xlabel='diameter (mm)',
        ylabel=ylabel, titl='snow density aloft',
        fname=fname, labels=label_list, logy=logy)


def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
