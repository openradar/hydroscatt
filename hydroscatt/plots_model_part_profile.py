#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
plots_model_part_profile
================================================

Auxiliary vertical profile plots of melting particles model

.. autosummary::
    :toctree: generated/

    plot_profile_d_ext
    plot_profile_d_int
    plot_profile_fmw
    plot_profile_mass
    plot_profile_snow_mass
    plot_profile_bulk_mass
    plot_profile_bulk_snow_mass
    plot_psd_at_altitude
    plot_psd_snow_at_altitude
    plot_psd_drops_at_altitude
    plot_fmw_at_altitude
    plot_vel_at_altitude
    plot_vel_snow_at_altitude
    plot_axis_ratio_at_altitude
    plot_nre_at_altitude

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
from warnings import warn
from copy import deepcopy

import numpy as np

from pytmatrix.psd import ExponentialPSD

from scattering_io import read_multiple_hydro_part_model, get_save_dir
from graph import plot_vertical_profile, plot_multiple_var
from part_descrip import compute_velocity_melting_hail, compute_reynolds_number
from part_descrip import compute_velocity_melting_snow, compute_equi_vol_diam
from atmos import compute_air_density, compute_dynamic_viscosity
from precip import psd_hail, compute_dsd_shed_water
from precip import compute_dsd_breakup_water

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to profile framework')

    # keyword arguments
    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='output data path')

    parser.add_argument(
        '--hydro_type', type=str,
        default='melting_hail',
        help='hydrometeor type. Default melting_hail')

    args = parser.parse_args()

    print(f'====== profile processing started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== profile processing finished: ")

    # list of initial hailstone diameter to plot
    d_init_vec = [4., 8., 10., 15., 25., 35.]
    alt_list = np.array([4000, 3000, 2000, 100, 0])  # masl
#    # alt_list = np.array([3300, 3330, 3350])  # masl
    ylabel = 'altitude [masl]'
    alt_label = 'masl'
    d_max = None

#    # list of initial snowflake diameter to plot
#    d_init_vec = [5., 10., 15., 20.]
#    alt_list = np.array([0, -100, -200, -300, -400, -500]) # m from iso0
#    d_max = 10.
#    ylabel = 'height with respect to iso0 [m]'
#    alt_label = 'm from iso0'

    df_model, temp_vec = read_multiple_hydro_part_model(
        args.path, args.hydro_type, d_max=d_max)

    if df_model is None:
        return

    savedir = get_save_dir(
        args.path, None, None, None, create_dir=True,
        with_subdirs=False)

    # Functions valid for all hydrometeors:
#    plot_profile_d_ext(
#        savedir, args.hydro_type, df_model, temp_vec, d_init_vec,
#        ylabel=ylabel)
#    plot_profile_d_int(
#        savedir, args.hydro_type, df_model, temp_vec, d_init_vec,
#        ylabel=ylabel)
#    plot_profile_fmw(
#        savedir, args.hydro_type, df_model, temp_vec, d_init_vec,
#        ylabel=ylabel)
#    plot_fmw_at_altitude(
#        savedir, args.hydro_type, df_model, temp_vec, alt_list,
#        alt_label=alt_label)
#    plot_axis_ratio_at_altitude(
#        savedir, args.hydro_type, df_model, temp_vec, alt_list,
#        alt_label=alt_label)

    if args.hydro_type == 'melting_hail':
#       plot_profile_mass(
#           savedir, args.hydro_type, df_model, temp_vec, d_init_vec)
#        plot_profile_bulk_mass(savedir, args.hydro_type, df_model, temp_vec)
#        plot_vel_at_altitude(
#            savedir, args.hydro_type, df_model, temp_vec, alt_list)
#        plot_nre_at_altitude(
#            savedir, args.hydro_type, df_model, temp_vec, alt_list)
#        plot_psd_at_altitude(
#            savedir, args.hydro_type, df_model, temp_vec, alt_list,
#            with_breakup=True)
        plot_psd_drops_at_altitude(
            savedir, args.hydro_type, df_model, temp_vec, alt_list,
            d_max_h=35., lamb_h=0.27, coeff_nh=800.)

#    if args.hydro_type == 'melting_snow':
#        plot_profile_snow_mass(
#            savedir, args.hydro_type, df_model, temp_vec, d_init_vec,
#            ylabel=ylabel)
#        plot_profile_bulk_snow_mass(
#            savedir, args.hydro_type, df_model, temp_vec, ylabel=ylabel)
#        plot_psd_snow_at_altitude(
#            savedir, args.hydro_type, df_model, temp_vec, alt_list)
#        plot_vel_snow_at_altitude(
#            savedir, args.hydro_type, df_model, temp_vec, alt_list)


def plot_profile_d_ext(savedir, hydro_type, df_model, temp_vec, d_init_vec,
                       lapse_rate=6.5, alt_iso0=4000.,
                       ylabel='altitude [masl]'):
    """
    Plots profile of external particle diameter

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    d_init_vec: array of float
        initial particle diameters to plot
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'

    Returns
    -------
    Nothing

    """
    d_ext_prof_list = []
    for d_init in d_init_vec:
        d_ext_prof = (
            df_model[np.isclose(df_model['d_init'], d_init)]['d_ext'].values)
        if d_ext_prof.size == 0:
            warn(f'no data for d_init {d_init}')
            continue
        d_ext_prof_list.append(d_ext_prof)

    if ylabel == 'altitude [masl]':
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False
    elif ylabel == 'temperature [deg C]':
        y_data = temp_vec
        invert_yaxis = True
    elif ylabel == 'height with respect to iso0 [m]':
        y_data = -temp_vec/lapse_rate*1000.
        invert_yaxis = False
    else:
        warn(f'unknown ylabel {ylabel}. Default altitude is going to be used')
        ylabel = 'altitude [masl]'
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False

    fname = f'{savedir}profile_{hydro_type}_alt-d_ext.png'
    plot_vertical_profile(
        d_ext_prof_list, y_data, xlabel='external particle diameter [mm]',
        ylabel=ylabel, titl=f'{hydro_type}', fname=fname,
        invert_yaxis=invert_yaxis)


def plot_profile_d_int(savedir, hydro_type, df_model, temp_vec, d_init_vec,
                       lapse_rate=6.5, alt_iso0=4000.,
                       ylabel='altitude [masl]'):
    """
    Plots profile of particle internal core diameter

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    d_init_vec: array of float
        initial particle diameters to plot
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'

    Returns
    -------
    Nothing

    """
    d_int_prof_list = []
    for d_init in d_init_vec:
        d_int_prof = (
            df_model[np.isclose(df_model['d_init'], d_init)]['d_int'].values)
        if d_int_prof.size == 0:
            warn(f'no data for d_init {d_init}')
            continue
        d_int_prof_list.append(d_int_prof)

    if ylabel == 'altitude [masl]':
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False
    elif ylabel == 'temperature [deg C]':
        y_data = temp_vec
        invert_yaxis = True
    elif ylabel == 'height with respect to iso0 [m]':
        y_data = -temp_vec/lapse_rate*1000.
        invert_yaxis = False
    else:
        warn(f'unknown ylabel {ylabel}. Default altitude is going to be used')
        ylabel = 'altitude [masl]'
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False

    fname = f'{savedir}profile_{hydro_type}_alt-d_int.png'
    plot_vertical_profile(
        d_int_prof_list, y_data, xlabel='internal particle diameter [mm]',
        ylabel=ylabel, titl=f'{hydro_type}', fname=fname,
        invert_yaxis=invert_yaxis)


def plot_profile_fmw(savedir, hydro_type, df_model, temp_vec, d_init_vec,
                     lapse_rate=6.5, alt_iso0=4000.,
                     ylabel='altitude [masl]'):
    """
    Plots profile of particle fraction mass water

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    d_init_vec: array of float
        initial particle diameters to plot
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'

    Returns
    -------
    Nothing

    """
    fmw_prof_list = []
    for d_init in d_init_vec:
        fmw_prof = (
            df_model[np.isclose(df_model['d_init'], d_init)]['fmw'].values)
        if fmw_prof.size == 0:
            warn(f'no data for d_init {d_init}')
            continue
        fmw_prof_list.append(fmw_prof)

    if ylabel == 'altitude [masl]':
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False
    elif ylabel == 'temperature [deg C]':
        y_data = temp_vec
        invert_yaxis = True
    elif ylabel == 'height with respect to iso0 [m]':
        y_data = -temp_vec/lapse_rate*1000.
        invert_yaxis = False
    else:
        warn(f'unknown ylabel {ylabel}. Default altitude is going to be used')
        ylabel = 'altitude [masl]'
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False

    fname = f'{savedir}profile_{hydro_type}_alt-fmw.png'
    plot_vertical_profile(
        fmw_prof_list, y_data, xlabel='mass water fraction',
        ylabel=ylabel, titl=f'{hydro_type}', fname=fname,
        labels=d_init_vec, invert_yaxis=invert_yaxis)


def plot_profile_mass(savedir, hydro_type, df_model, temp_vec, d_init_vec,
                      lapse_rate=6.5, alt_iso0=4000.,
                      ylabel='altitude [masl]'):
    """
    Plots profile of hailstone mass

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    d_init_vec: array of float
        initial hailstone diameters to plot
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'

    Returns
    -------
    Nothing

    """
    for d_init in d_init_vec:
        hail_mass = (
            df_model[np.isclose(df_model['d_init'], d_init)][
                'hail_mass'].values)
        shed_water_mass = (
            df_model[np.isclose(df_model['d_init'], d_init)][
                'shed_water_mass'].values)
        fmw = df_model[np.isclose(df_model['d_init'], d_init)]['fmw'].values

        water_mass = hail_mass*fmw
        ice_mass = hail_mass - water_mass

        if ylabel == 'altitude [masl]':
            y_data = alt_iso0-temp_vec/lapse_rate*1000.
            invert_yaxis = False
        elif ylabel == 'temperature [deg C]':
            y_data = temp_vec
            invert_yaxis = True
        elif ylabel == 'height with respect to iso0 [m]':
            y_data = -temp_vec/lapse_rate*1000.
            invert_yaxis = False
        else:
            warn(f'unknown ylabel {ylabel}. '
                 f'Default altitude is going to be used')
            ylabel = 'altitude [masl]'
            y_data = alt_iso0-temp_vec/lapse_rate*1000.
            invert_yaxis = False

        fname = (
            f'{savedir}profile_{hydro_type}'
            f'_d{int(d_init):02d}-mass.png')
        plot_vertical_profile(
            [ice_mass, water_mass, shed_water_mass],
            y_data, xlabel='mass (g)', ylabel='altitude (masl)',
            titl=f'{hydro_type} d_init {d_init} mm', fname=fname,
            labels=['ice', 'melted water', 'shed water'],
            invert_yaxis=invert_yaxis)


def plot_profile_snow_mass(savedir, hydro_type, df_model, temp_vec,
                           d_init_vec, lapse_rate=6.5, alt_iso0=4000.,
                           ylabel='altitude [masl]'):
    """
    Plots profile of snow mass

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    d_init_vec: array of float
        initial hailstone diameters to plot
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'

    Returns
    -------
    Nothing

    """
    for d_init in d_init_vec:
        snow_mass = (
            df_model[np.isclose(df_model['d_init'], d_init)][
                'mass'].values)
        fmw = df_model[np.isclose(df_model['d_init'], d_init)]['fmw'].values

        water_mass = snow_mass*fmw
        ice_mass = snow_mass - water_mass

        if ylabel == 'altitude [masl]':
            y_data = alt_iso0-temp_vec/lapse_rate*1000.
            invert_yaxis = False
        elif ylabel == 'temperature [deg C]':
            y_data = temp_vec
            invert_yaxis = True
        elif ylabel == 'height with respect to iso0 [m]':
            y_data = -temp_vec/lapse_rate*1000.
            invert_yaxis = False
        else:
            warn(f'unknown ylabel {ylabel}.'
                 f' Default altitude is going to be used')
            ylabel = 'altitude [masl]'
            y_data = alt_iso0-temp_vec/lapse_rate*1000.
            invert_yaxis = False

        fname = (
            f'{savedir}profile_{hydro_type}'
            f'_d{int(d_init):02d}-mass.png')
        plot_vertical_profile(
            [ice_mass, water_mass],
            y_data, xlabel='mass (g)', ylabel='altitude (masl)',
            titl=f'{hydro_type} d_init {d_init} mm', fname=fname,
            labels=['ice', 'melted water'],
            invert_yaxis=invert_yaxis)


def plot_profile_bulk_mass(savedir, hydro_type, df_model, temp_vec,
                           lapse_rate=6.5, alt_iso0=4000.,
                           ylabel='altitude [masl]'):
    """
    Plots profile of bulk mass as a function of altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'

    Returns
    -------
    Nothing

    """
    ice_mass = np.empty(temp_vec.size)
    water_mass = np.empty(temp_vec.size)
    shed_water_mass = np.empty(temp_vec.size)
    break_up_mass = np.empty(temp_vec.size)

    for ind, temp in enumerate(temp_vec):
        df_aux = df_model[np.isclose(df_model['temp'], temp)].copy()

        df_aux.sort_values('d_init', inplace=True)
        d_init = df_aux['d_init'].values
        d_ext = df_aux['d_ext'].values
        alpha = df_aux['alpha'].values
        hail_mass_aux = df_aux['hail_mass'].values
        fmw_aux = df_aux['fmw'].values
        shed_water_mass_aux = df_aux['shed_water_mass'].values
        prob_break = df_aux['prob_break'].values
        vel_hail = df_aux['vel'].values
        water_mass_aux = hail_mass_aux*fmw_aux
        ice_mass_aux = hail_mass_aux - water_mass_aux

        bin_left = np.append([0], d_init[:-1])
        bin_right = deepcopy(d_init)
        delta_d = bin_right - bin_left

        if np.isclose(temp, 0.):
            vel_hail0 = deepcopy(vel_hail)

        psd_hail_no_break = psd_hail(d_init)

        # conservation of the flux
        psd_hail_no_break *= (vel_hail0/vel_hail)

        psd_break_vals = psd_hail_no_break*prob_break
        psd_hail_vals = psd_hail_no_break-psd_break_vals

        break_up_mass[ind] = np.sum(water_mass_aux*psd_break_vals*delta_d)
        water_mass[ind] = np.sum(water_mass_aux*psd_hail_vals*delta_d)
        ice_mass[ind] = np.sum(ice_mass_aux*psd_hail_vals*delta_d)
        shed_water_mass[ind] = np.sum(
            shed_water_mass_aux*psd_hail_no_break*delta_d)

    total_mass = ice_mass+water_mass+shed_water_mass+break_up_mass

    if ylabel == 'altitude [masl]':
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False
    elif ylabel == 'temperature [deg C]':
        y_data = temp_vec
        invert_yaxis = True
    elif ylabel == 'height with respect to iso0 [m]':
        y_data = -temp_vec/lapse_rate*1000.
        invert_yaxis = False
    else:
        warn(f'unknown ylabel {ylabel}. Default altitude is going to be used')
        ylabel = 'altitude [masl]'
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False

    fname = (
        f'{savedir}profile_{hydro_type}'
        f'_alt-bulk-mass.png')
    plot_vertical_profile(
        [ice_mass, water_mass, shed_water_mass, break_up_mass, total_mass],
        y_data, xlabel='mass content (g/m3)', ylabel=ylabel,
        titl=f'{hydro_type}', fname=fname,
        labels=['ice', 'melted water', 'shed water', 'breakup water',
                'total_mass'],
        invert_yaxis=invert_yaxis)


def plot_profile_bulk_snow_mass(savedir, hydro_type, df_model, temp_vec,
                                lapse_rate=6.5, alt_iso0=4000.,
                                ylabel='altitude [masl]', nw=8000.,
                                rr=5.):
    """
    Plots profile of bulk mass as a function of altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec: array of floats
        temperature values of the profile
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    ylabel : str
        the label defining the y-axis. Can be 'altitude [masl]',
        'temperature [deg C]' or 'height with respect to iso0 [m]'
    nw : float
        N0 value of an exponentinal distribution mm^-1 m^-3
    rr : float
        liquid equivalent rainfall rate (mm/h)

    Returns
    -------
    Nothing

    """
    lamb = 4.1*np.power(rr, -0.21)

    ice_mass = np.empty(temp_vec.size)
    water_mass = np.empty(temp_vec.size)

    for ind, temp in enumerate(temp_vec):
        df_aux = df_model[np.isclose(df_model['temp'], temp)].copy()

        df_aux.sort_values('d_init', inplace=True)
        d_init = df_aux['d_init'].values
        snow_mass_aux = df_aux['mass'].values
        fmw_aux = df_aux['fmw'].values
        vel = df_aux['vel'].values
        if np.isclose(temp, 0.):
            vel0 = deepcopy(vel)
        
        water_mass_aux = snow_mass_aux*fmw_aux
        ice_mass_aux = snow_mass_aux - water_mass_aux
        

        d_rd = compute_equi_vol_diam(snow_mass_aux)
        psd_rain = ExponentialPSD(
            N0=nw, Lambda=lamb, D_max=d_rd[-1])
        psd_vals_rain = psd_rain(d_rd)
        delta_d_rain = d_rd - np.append(0, d_rd[:-1])

        bin_left = np.append([0], d_init[:-1])
        bin_right = deepcopy(d_init)
        delta_d_snow = bin_right - bin_left

        psd_vals_snow = psd_vals_rain*delta_d_rain/delta_d_snow
        
        # conservation of the flux
        psd_vals_snow *= (vel0/vel)

        water_mass[ind] = np.sum(water_mass_aux*psd_vals_snow*delta_d_snow)
        ice_mass[ind] = np.sum(ice_mass_aux*psd_vals_snow*delta_d_snow)

    total_mass = ice_mass+water_mass

    if ylabel == 'altitude [masl]':
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False
    elif ylabel == 'temperature [deg C]':
        y_data = temp_vec
        invert_yaxis = True
    elif ylabel == 'height with respect to iso0 [m]':
        y_data = -temp_vec/lapse_rate*1000.
        invert_yaxis = False
    else:
        warn(f'unknown ylabel {ylabel}. Default altitude is going to be used')
        ylabel = 'altitude [masl]'
        y_data = alt_iso0-temp_vec/lapse_rate*1000.
        invert_yaxis = False

    fname = (
        f'{savedir}profile_{hydro_type}'
        f'_alt-bulk-mass.png')
    plot_vertical_profile(
        [ice_mass, water_mass, total_mass],
        y_data, xlabel='mass content (g/m3)', ylabel=ylabel,
        titl=f'{hydro_type}', fname=fname,
        labels=['ice', 'melted water', 'total_mass'],
        invert_yaxis=invert_yaxis)


def plot_psd_at_altitude(savedir, hydro_type, df_model, temp_vec, alt_list,
                         d_max_h=24., lamb_h=0.42, coeff_nh=400.,
                         lapse_rate=6.5, alt_iso0=4000., with_breakup=True):
    """
    Plots hail PSD at given altitudes

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot (km)
    dmax_h, lamb_h, coeff_nh : float
        parameters that define the hail PSD
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    with_breakup : bool
        if True breakup of large drops is considered

    Returns
    -------
    Nothing

    """
    for alt in alt_list:
        temp = (alt_iso0-alt)*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        alt_ref = alt_iso0-temp_ref/lapse_rate*1000.

        df_aux = df_model[
            np.isclose(df_model['temp'], temp_ref)].copy()

        df_aux.sort_values('d_init', inplace=True)
        d_init = df_aux['d_init'].values
        d_int = df_aux['d_int'].values
        d_ext = df_aux['d_ext'].values
        prob_break = df_aux['prob_break'].values
        vel_hail = df_aux['vel'].values

        bin_left = np.append([0], d_init[:-1])
        bin_right = deepcopy(d_init)

        if np.isclose(alt, alt_iso0):
            vel_hail0 = deepcopy(vel_hail)

        psd_hail_no_break = psd_hail(
            d_init, lamb_h=lamb_h, coeff_nh=coeff_nh, d_max_h=d_max_h)

        # conservation of the flux
        psd_hail_no_break *= (vel_hail0/vel_hail)

        if with_breakup:
            psd_break_vals = psd_hail_no_break*prob_break
            label = 'with_breakup'
        else:
            psd_break_vals = np.zeros(psd_hail_no_break.size)
            label = 'no_breakup'

        psd_hail_vals_aux = psd_hail_no_break-psd_break_vals

        psd_hail_vals_aux[psd_hail_vals_aux == 0.] = np.nan
        psd_ice_core_aux = deepcopy(psd_hail_vals_aux)
        psd_ice_core_aux[d_int <= 0.001] = np.nan

        bin_left_vec = np.append([0], d_init[:-1])
        bin_right_vec = deepcopy(d_init)

        psd_ice_core = np.zeros(d_init.size)
        psd_hail_vals = np.zeros(d_init.size)
        for ind_bin, (bin_left, bin_right) in enumerate(zip(
                bin_left_vec, bin_right_vec)):
            ind_close = np.where(
                (d_ext > bin_left) & (d_ext <= bin_right))[0]
            if ind_close.size > 0:
                psd_hail_vals[ind_bin] = np.sum(
                    psd_hail_vals_aux[ind_close])

            ind_close = np.where(
                (d_int > bin_left) & (d_int <= bin_right))[0]
            if ind_close.size > 0:
                psd_ice_core[ind_bin] = np.sum(
                    psd_ice_core_aux[ind_close])
        psd_hail_vals[psd_hail_vals == 0] = np.nan
        psd_ice_core[psd_ice_core == 0] = np.nan

        fname = (
            f'{savedir}profile_{hydro_type}_{label}_alt{alt_ref:.0f}-psd.png')
        plot_multiple_var(
            [d_init, d_init], [psd_ice_core, psd_hail_vals],
            xlabel='diameter (mm)', ylabel='N(D) (m-3 mm-1)',
            titl=f'{hydro_type} {label} PSD at {alt_ref:.0f} masl',
            fname=fname, labels=['ice core', 'hailstone'], logy=True)


def plot_psd_snow_at_altitude(savedir, hydro_type, df_model, temp_vec,
                              alt_list, nw=9370., lamb=1.18, lapse_rate=6.5):
    """
    Plots snow PSD at given altitudes

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes with respect to iso0 to plot
    nw, lamb : float
        parameters that define the snow PSD
    lapse_rate : float
        temperature lapse rate deg/km

    Returns
    -------
    Nothing

    """
    for alt in alt_list:
        temp = -alt*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        alt_ref = -temp_ref/lapse_rate*1000.

        df_aux = df_model[
            np.isclose(df_model['temp'], temp_ref)].copy()

        df_aux.sort_values('d_init', inplace=True)
        d_init = df_aux['d_init'].values
        d_ext = df_aux['d_ext'].values
        vel = df_aux['vel'].values

        bin_left = np.append([0], d_init[:-1])
        bin_right = deepcopy(d_init)
        
        if np.isclose(alt, 0):
            vel0 = deepcopy(vel)

        psd_snow = ExponentialPSD(N0=nw, Lambda=lamb, D_max=d_init[-1])
        psd_vals_snow_aux = psd_snow(d_init)
        
        # conservation of the flux
        psd_vals_snow_aux *= (vel0/vel)

        bin_left_vec = np.append([0], d_init[:-1])
        bin_right_vec = deepcopy(d_init)

        psd_vals_snow = np.zeros(d_init.size)
        for ind_bin, (bin_left, bin_right) in enumerate(zip(
                bin_left_vec, bin_right_vec)):
            ind_close = np.where(
                (d_ext > bin_left) & (d_ext <= bin_right))[0]
            if ind_close.size > 0:
                psd_vals_snow[ind_bin] = np.sum(
                    psd_vals_snow_aux[ind_close])

        psd_vals_snow[psd_vals_snow == 0] = np.nan

        fname = (
            f'{savedir}profile_{hydro_type}_alt{alt_ref:.0f}-psd.png')
        plot_multiple_var(
            [d_init], [psd_vals_snow],
            xlabel='diameter (mm)', ylabel='N(D) (m-3 mm-1)',
            titl=f'{hydro_type} PSD at {alt_ref:.0f} m respect to iso0',
            fname=fname, logy=True)


def plot_psd_drops_at_altitude(savedir, hydro_type, df_model, temp_vec,
                               alt_list, d_max_h=24., lamb_h=0.42,
                               coeff_nh=400., lapse_rate=6.5, alt_iso0=4000.):
    """
    Plots hail PSD at given altitudes

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot (masl)
    dmax_h, lamb_h, coeff_nh : float
        parameters that define the hail PSD
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)

    Returns
    -------
    Nothing

    """
    # parameters for rainfall created from shed water
    diam_rd_min = 0.1
    diam_rd_max = 4.5
    step = 0.1
    diam_rd = np.arange(diam_rd_min, diam_rd_max+step, step)

    for alt in alt_list:
        temp = (alt_iso0-alt)*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        alt_ref = alt_iso0-temp_ref/lapse_rate*1000.

        df_aux = df_model[
            np.isclose(df_model['temp'], temp_ref)].copy()

        df_aux.sort_values('d_init', inplace=True)
        d_init = df_aux['d_init'].values
        d_ext = df_aux['d_ext'].values
        fmw = df_aux['fmw'].values
        shed_water_mass = df_aux['shed_water_mass'].values
        water_mass = df_aux['hail_mass'].values*fmw
        prob_break = df_aux['prob_break'].values
        vel_hail = df_aux['vel'].values

        if np.isclose(alt, alt_iso0):
            vel_hail0 = deepcopy(vel_hail)

        bin_left = np.append([0], d_init[:-1])
        bin_right = deepcopy(d_init)

        psd_hail_no_break = psd_hail(
            d_init, lamb_h=lamb_h, coeff_nh=coeff_nh, d_max_h=d_max_h)
        # conservation of the flux
        psd_hail_no_break *= (vel_hail0/vel_hail)

        psd_break_vals = psd_hail_no_break*prob_break
        psd_hail_vals_aux = psd_hail_no_break-psd_break_vals

        psd_hail_vals_aux[psd_hail_vals_aux == 0.] = np.nan
        psd_melted_aux = deepcopy(psd_hail_vals_aux)
        psd_melted_aux[fmw < 1.] = np.nan

        psd_shed = np.zeros(diam_rd.size)
        if np.sum(shed_water_mass) > 0:
            psd_shed_func = compute_dsd_shed_water(
                d_init, psd_hail_no_break, diam_rd, shed_water_mass)
            psd_shed = psd_shed_func(diam_rd)

        d_breakup = d_ext[np.isclose(fmw, 1.)]
        psd_break = np.zeros(d_breakup.size)
        if np.sum(prob_break) > 0:
            psd_break_func = compute_dsd_breakup_water(
                d_init, psd_break_vals, d_breakup, water_mass)
            psd_break = psd_break_func(d_breakup)

        bin_left_vec = np.append([0], d_init[:-1])
        bin_right_vec = deepcopy(d_init)

        psd_melted = np.zeros(d_init.size)
        for ind_bin, (bin_left, bin_right) in enumerate(zip(
                bin_left_vec, bin_right_vec)):
            ind_close = np.where(
                (d_ext > bin_left) & (d_ext <= bin_right))[0]
            if ind_close.size > 0:
                psd_melted[ind_bin] = np.sum(psd_melted_aux[ind_close])
        psd_break[psd_break == 0] = np.nan
        psd_shed[psd_shed == 0] = np.nan
        psd_melted[psd_melted == 0] = np.nan

        fname = (
            f'{savedir}profile_{hydro_type}_alt{alt_ref:.0f}-dsd.png')
        plot_multiple_var(
            [diam_rd, d_breakup, d_init], [psd_shed, psd_break, psd_melted],
            xlabel='diameter (mm)', ylabel='N(D) (m3 mm-1)',
            titl=f'{hydro_type} DSD at {alt_ref:.0f} masl', fname=fname,
            labels=['shed', 'breakup', 'melted'], logy=True)


def plot_fmw_at_altitude(savedir, hydro_type, df_model, temp_vec, alt_list,
                         lapse_rate=6.5, alt_iso0=4000., alt_label='masl'):
    """
    Plots fmw at given altitudes

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot (masl or m with respect to iso0)
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    alt_label : str
        altitude label can be 'masl' or 'm from iso0'

    Returns
    -------
    Nothing

    """
    fmw_temp_list = []
    d_ext_temp_list = []
    label_list = []
    for alt in alt_list:
        if alt_label == 'masl':
            temp = (alt_iso0-alt)*lapse_rate/1000.
        else:
            temp = -alt*lapse_rate/1000.

        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        if alt_label == 'masl':
            alt_ref = alt_iso0-temp_ref/lapse_rate*1000.
        else:
            alt_ref = -temp_ref/lapse_rate*1000.

        df_aux = df_model[np.isclose(df_model['temp'], temp_ref)].copy()
        df_aux.sort_values('d_init', inplace=True)
        d_ext = df_aux['d_ext'].values
        fmw = df_aux['fmw'].values
        if fmw.size == 0:
            warn(f'no data for temp {temp}')
            continue
        fmw_temp_list.append(fmw)
        d_ext_temp_list.append(d_ext)

        label_list.append(f'{alt_ref:.0f} {alt_label} - {temp_ref:.1f} degC')

    fname = f'{savedir}profile_{hydro_type}_d_ext-fmw.png'
    plot_multiple_var(
        d_ext_temp_list, fmw_temp_list,
        xlabel='external hailstone diameter (mm)',
        ylabel='mass water fraction', titl=f'{hydro_type}', fname=fname,
        labels=label_list)


def plot_vel_at_altitude(savedir, hydro_type, df_model, temp_vec, alt_list,
                         lapse_rate=6.5, alt_iso0=4000.):
    """
    Plots hailstone terminal velocity at a given altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot (masl)
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)

    Returns
    -------
    Nothing

    """
    for alt in alt_list:
        temp = (alt_iso0-alt)*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        alt_ref = alt_iso0-temp_ref/lapse_rate*1000.

        df_aux = df_model[np.isclose(df_model['temp'], temp_ref)].copy()
        df_aux.sort_values('d_ext', inplace=True)
        d_ext = df_aux['d_ext'].values
        fmw = df_aux['fmw'].values
        alpha = df_aux['alpha'].values
        mass = df_aux['hail_mass'].values
        if fmw.size == 0:
            warn(f'no data for temp {temp_ref}')
            continue
        vel, _ = compute_velocity_melting_hail(
            d_ext, mass, fmw, alpha, dens_air=compute_air_density(alt_ref),
            visc=compute_dynamic_viscosity(temp_ref))

        fname = (
            f'{savedir}profile_{hydro_type}_{int(temp_ref*100):04d}'
            f'_d_ext-vel.png')
        titl = (
            f'{hydro_type} velocity at {temp_ref:.2f} degC - {alt_ref} masl')
        plot_multiple_var(
            [d_ext], [vel], xlabel='external hailstone diameter (mm)',
            ylabel='hailstone terminal velocity (m/s)',
            titl=titl, fname=fname)


def plot_vel_snow_at_altitude(savedir, hydro_type, df_model, temp_vec,
                              alt_list, lapse_rate=6.5):
    """
    Plots hailstone terminal velocity at a given altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot (masl)
    lapse_rate : float
        temperature lapse rate deg/km

    Returns
    -------
    Nothing

    """
    for alt in alt_list:
        temp = -alt*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        alt_ref = -temp_ref/lapse_rate*1000.

        df_aux = df_model[np.isclose(df_model['temp'], temp_ref)].copy()
        df_aux.sort_values('d_ext', inplace=True)
        d_ext = df_aux['d_ext'].values
        fmw = df_aux['fmw'].values
        mass = df_aux['mass'].values
        if fmw.size == 0:
            warn(f'no data for temp {temp_ref}')
            continue
        vel = compute_velocity_melting_snow(d_ext, mass, fmw)

        fname = (
            f'{savedir}profile_{hydro_type}_{int(temp_ref*100):04d}'
            f'_d_ext-vel.png')
        titl = (
            f'{hydro_type} velocity at {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} m from iso0')
        plot_multiple_var(
            [d_ext], [vel], xlabel='snowflake diameter (mm)',
            ylabel='snowflake terminal velocity (m/s)',
            titl=titl, fname=fname)


def plot_axis_ratio_at_altitude(savedir, hydro_type, df_model, temp_vec,
                                alt_list, lapse_rate=6.5, alt_iso0=4000.,
                                alt_label='masl'):
    """
    Plots hailstone terminal velocity at a given altitude

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot (masl)
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude
    alt_label : str
        altitude reference label

    Returns
    -------
    Nothing

    """
    for alt in alt_list:
        if alt_label == 'masl':
            temp = (alt_iso0-alt)*lapse_rate/1000.
        else:
            temp = -alt*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        if alt_label == 'masl':
            alt_ref = alt_iso0-temp_ref/lapse_rate*1000.
        else:
            alt_ref = -temp_ref/lapse_rate*1000.

        df_aux = df_model[np.isclose(df_model['temp'], temp_ref)].copy()
        df_aux.sort_values('d_ext', inplace=True)
        d_ext = df_aux['d_ext'].values
        d_int = df_aux['d_int'].values
        ar_ext = df_aux['ar_ext'].values
        ar_int = df_aux['ar_int'].values
        if d_ext.size == 0:
            warn(f'no data for temp {temp_ref}')
            continue

        fname = (
            f'{savedir}profile_{hydro_type}_{int(temp_ref*100):04d}'
            f'_d_ext-ar.png')
        titl = (
            f'{hydro_type} axis ratio at {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} {alt_label}')
        plot_multiple_var(
            [d_ext], [ar_ext], xlabel='external particle diameter (mm)',
            ylabel='axis ratio', titl=titl, fname=fname)

        fname = (
            f'{savedir}profile_{hydro_type}_{int(temp_ref*100):04d}'
            f'_d_int-ar.png')
        titl = (
            f'{hydro_type} axis ratio at {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} {alt_label}')
        plot_multiple_var(
            [d_int], [ar_int], xlabel='internal core diameter (mm)',
            ylabel='axis ratio', titl=titl, fname=fname)


def plot_nre_at_altitude(savedir, hydro_type, df_model, temp_vec, alt_list,
                         lapse_rate=6.5, alt_iso0=4000.):
    """
    Plots hailstone Reynolds number at different altitudes

    Parameters
    ----------
    savedir : str
        directory where to save the data
    hydro_type: str
        hydrometeor type
    df_model: pandas dataframe
        dataframe containing the data
    temp_vec : array of float
        temperature values of the profile
    alt_list: array of float
        list of altitudes to plot
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)

    Returns
    -------
    Nothing

    """
    for alt in alt_list:
        temp = (alt_iso0-alt)*lapse_rate/1000.
        ind_temp = np.argmin(np.abs(temp_vec-temp))
        temp_ref = temp_vec[ind_temp]
        alt_ref = alt_iso0-temp_ref/lapse_rate*1000.

        df_aux = df_model[np.isclose(df_model['temp'], temp_ref)].copy()
        df_aux.sort_values('d_ext', inplace=True)
        d_ext = df_aux['d_ext'].values
        fmw = df_aux['fmw'].values
        mass = df_aux['hail_mass'].values
        if fmw.size == 0:
            warn(f'no data for temp {temp_ref}')
            continue
        nre = compute_reynolds_number(
            mass, dens_air=compute_air_density(alt_ref),
            visc=compute_dynamic_viscosity(temp_ref))

        fname = (
            f'{savedir}profile_{hydro_type}_{int(temp_ref*100):04d}'
            f'_d_ext-nre.png')
        titl = (
            f'{hydro_type} Reynolds number at {temp_ref:.2f} degC '
            '- {alt_ref} masl')
        plot_multiple_var(
            [d_ext], [nre], xlabel='external hailstone diameter (mm)',
            ylabel='Reynolds number', titl=titl, fname=fname)


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
