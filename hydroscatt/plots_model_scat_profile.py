#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
plots_model_scat_profile
================================================

Auxiliary vertical profile plots of scattering model

.. autosummary::
    :toctree: generated/

    plot_diel_constant_at_altitude

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
from warnings import warn

import numpy as np

from scattering_io import read_multiple_hydro_scatt_model, get_save_dir
from graph import plot_multiple_var


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

    parser.add_argument(
        '--band', type=str,
        default='C',
        help='Frequency band. Default C')

    args = parser.parse_args()

    print(f'====== profile processing started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== profile processing finished: ")

    # list of hail altitudes to plot
    # alt_list = np.array([0, 100, 2000, 3000, 4000])  # masl
    alt_list = np.array([3300, 3330, 3350])  # masl
    alt_label = 'masl'

    # list of snow altitudes to plot
    # alt_list = np.array([0, -100, -200, -300, -400, -500]) # m from iso0
    alt_label = 'm from iso0'

    df_model, temp_vec = read_multiple_hydro_scatt_model(
        args.path, args.hydro_type, band=args.band)

    if df_model is None:
        return

    savedir = get_save_dir(
        args.path, None, None, None, create_dir=True,
        with_subdirs=False)

    plot_diel_constant_at_altitude(
         savedir, args.hydro_type, df_model, temp_vec, alt_list,
         band=args.band, alt_label=alt_label)


def plot_diel_constant_at_altitude(savedir, hydro_type, df_model, temp_vec,
                                   alt_list, band='C', lapse_rate=6.5,
                                   alt_iso0=4000., alt_label='masl'):
    """
    Plots dielectric constant at a given altitude

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
    band : str
        frequency band at which the dielectric constant was computed
    lapse_rate : float
        temperature lapse rate deg/km
    alt_iso0 : float
        iso0 altitude (masl)
    alt_label : str
        the label defining the altitude. Can be 'masl' or 'm from iso0'

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
        d_ext = df_aux['d_ext'].values*10.
        d_int = df_aux['d_int'].values*10.
        err_ext = df_aux['err_ext'].values
        eri_ext = df_aux['eri_ext'].values
        err_int = df_aux['err_int'].values
        eri_int = df_aux['eri_int'].values

        if d_ext.size == 0:
            warn(f'no data for temp {temp_ref}')
            continue

        d_int[d_int == 0.001] = np.nan

        fname = (
            f'{savedir}profile_{hydro_type}_{band}_{int(temp_ref*100):04d}'
            f'_d_int-diel_re.png')
        titl = (
            f'{hydro_type} Re(dielectric constant)'
            f'\n{band}-band {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} {alt_label}')
        plot_multiple_var(
            [d_int], [err_int],
            xlabel='internal diameter (mm)', ylabel='dielectric constant',
            titl=titl, fname=fname)

        fname = (
            f'{savedir}profile_{hydro_type}_{band}_{int(temp_ref*100):04d}'
            '_d_int-diel_im.png')
        titl = (
            f'{hydro_type} Im(dielectric constant)\n'
            f'{band}-band {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} {alt_label}')
        plot_multiple_var(
            [d_int], [eri_int],
            xlabel='internal diameter (mm)', ylabel='dielectric constant',
            titl=titl, fname=fname)

        fname = (
            f'{savedir}profile_{hydro_type}_{band}_{int(temp_ref*100):04d}'
            f'_d_ext-diel_re.png')
        titl = (
            f'{hydro_type} Re(dielectric constant)\n'
            f'{band}-band {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} {alt_label}')
        plot_multiple_var(
            [d_ext], [err_ext],
            xlabel='external diameter (mm)', ylabel='dielectric constant',
            titl=titl, fname=fname)

        fname = (
            f'{savedir}profile_{hydro_type}_{band}_{int(temp_ref*100):04d}'
            f'_d_ext-diel_im.png')
        titl = (
            f'{hydro_type} Im(dielectric constant)'
            f'\n{band}-band {temp_ref:.2f} degC'
            f' - {alt_ref:.0f} {alt_label}')
        plot_multiple_var(
            [d_ext], [eri_ext],
            xlabel='external diameter (mm)', ylabel='dielectric constant',
            titl=titl, fname=fname)


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
