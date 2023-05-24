#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
scattering_snow
================================================

modelling of scattering properties of snow

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
import os
from warnings import warn
from copy import deepcopy

import numpy as np
import pandas as pd

from pytmatrix.psd import ExponentialPSD

from scattering_io import read_melting_hydro_part_model
from scattering_io import read_scatt_double_layer, read_wavelength_file

from part_descrip import compute_equi_vol_diam

from precip import compute_elwc, compute_equi_rainfall_rate

from refractivity import wavelength_to_band

from scattering import compute_angular_moments, compute_scattering_canting_psd
from scattering import compute_angular_moments_analytical

from graph import plot_psd_scatt_profile

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to snow scattering simulations framework')

    # keyword arguments
    parser.add_argument(
        '--input_path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/profile_snow/',
        help='input data path')

    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/profile_snow/',
        help='data path')

    parser.add_argument(
        '--band', type=str,
        default='C',
        help='frequency band. Default C')

    parser.add_argument(
        '--ele', type=float,
        default=0.,
        help='elevation angle. Default 0 deg')

    parser.add_argument(
        '--rr', type=float,
        default=5.,
        help='equivalent rainfall rate. Default 5 mm/h')
        
    parser.add_argument(
        '--d_max', type=float,
        default=None,
        help='Maximum snowflake size. If None it will be that of the input file. Default None')

    parser.add_argument(
        '--analytical_cant_angl', type=int,
        default=1,
        help='If 1 the canting angle will be computed analytically. Default 1')

    args = parser.parse_args()

    print(f'====== snow scattering simulation started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== snow scattering simulation finished: ")

    if args.analytical_cant_angl == 1 and args.ele != 0:
        warn('The analytical computation of canting angle is only valid for '
             'elevation angles close to 0 deg')
        return

    hydro_type = 'melting_snow'

    # Exponential distribution according to Gunn and Marshall
    # Ns(Dw) = N0 exp(-lamb*Dw)
    # N0 = 3.8e3*R^-0.87
    # lamb = 2.55*R^-0.48
    hydro_label = f'rreq{args.rr:.1f}_melting_snow'
    nw = 8000.
    lamb = 4.1*np.power(args.rr, -0.21)

    # parameters
    psd_var_list = [
        'refl_h', 'refl_v', 'ldr_h', 'ldr_v', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A_h', 'A_v', 'Adp']

    psd_x_var_list = [
        'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']
    psd_y_var_list = ['temp']

    print(hydro_type)
    print(hydro_label)
    print('n0', nw)
    print('lamb', lamb)
    print('band', args.band)
    print('elevation angle', args.ele)

    flist_model = glob.glob(
        f'{args.input_path}sp_{hydro_type}_*_model_part.csv')
    flist_scatt = glob.glob(
        f'{args.input_path}sp_{hydro_type}_{args.band}_*_tmat.out')
    if not flist_model or not flist_scatt:
        if not flist_model:
            warn(f'no model file at '
                 f'{args.input_path}sp_{hydro_type}_*_model_part.csv')
        if not flist_scatt:
            warn(f'no scattering file at '
                 f'{args.input_path}sp_{hydro_type}_{args.band}'
                 f'_*_tmat.out')
        return
    if len(flist_model) != len(flist_scatt):
        warn(f'Number of model files {len(flist_model)} different from '
             f'number of scattering files {len(flist_scatt)}')
        return

    freq_file = (
        f'{args.input_path}sp_{hydro_type}_{args.band}_freq.inp')
    wavelength = read_wavelength_file(freq_file)
    band = wavelength_to_band(wavelength)

    ntemp = len(flist_model)
    psd_snow_dict = {
        'temp': np.zeros(ntemp),
        'alt': np.zeros(ntemp),
        'lwc': np.zeros(ntemp),
        'rr': np.zeros(ntemp)
    }
    for var in psd_var_list:
        psd_snow_dict.update({var: np.zeros(ntemp)})

    # output file name
    fname = (
        f'{args.path}psd_profile_{hydro_label}_{band}'
        f'_ele{int(args.ele*100.):05d}_scattering.csv')
    for ind_temp, (model_file, scatt_file) in enumerate(zip(
            flist_model, flist_scatt)):
        df_model, temp = read_melting_hydro_part_model(
            model_file, d_max=args.d_max)
        n_sf = df_model.shape[0]
        df_scatt = read_scatt_double_layer(scatt_file, nrows=n_sf)

        psd_snow_dict['temp'][ind_temp] = temp
        print('temp', temp)
        print('computing angular moments ...')
        canting_angle_snow = 40.+df_model['fmw'].values*(10-40)
        if args.analytical_cant_angl:
            ang_moments_dict = compute_angular_moments_analytical(
                canting_angle_snow)
        else:
            ang_moments_dict = compute_angular_moments(
                canting_angle_snow, ele=args.ele)

        print('computing PSD scattering parameters of melting snowflakes ...')

        # equivalent raindrop diameter
        d_rd = compute_equi_vol_diam(df_model['mass'].values)
        psd_rain = ExponentialPSD(
            N0=nw, Lambda=lamb, D_max=d_rd[-1])
        psd_vals_rain = psd_rain(d_rd)
        delta_d_rain = d_rd - np.append(0, d_rd[:-1])

        delta_d_snow = (
            df_model['d_init'].values
            - np.append(0, df_model['d_init'].values[:-1]))

        vel_snow = df_model['vel'].values

        if np.isclose(temp, 0):
            vel_snow0 = deepcopy(vel_snow)

        psd_vals_snow = psd_vals_rain*delta_d_rain/delta_d_snow

        # conservation of the flux
        psd_vals_snow *= (vel_snow0/vel_snow)

        psd_snow_dict['lwc'][ind_temp] = compute_elwc(
            delta_d_snow, df_model['mass'].values, psd_vals_snow)
        psd_snow_dict['rr'][ind_temp] = compute_equi_rainfall_rate(
            delta_d_snow, df_model['mass'].values, psd_vals_snow,
            vel_snow)

        snow_psd_dict = compute_scattering_canting_psd(
            wavelength, df_scatt['fv180'].values, df_scatt['fh180'].values,
            df_scatt['fv0'].values, df_scatt['fh0'].values, ang_moments_dict,
            delta_d_snow, psd_vals_snow, var_list=psd_var_list)

        for var in psd_var_list:
            psd_snow_dict[var][ind_temp] = snow_psd_dict[var]

        # save DSD scattering results
        psd_snow_dict_aux = {}
        for key, value in psd_snow_dict.items():
            psd_snow_dict_aux.update({key: [value[ind_temp]]})
        df_snow_psd_aux = pd.DataFrame.from_dict(psd_snow_dict_aux)
        df_snow_psd_aux.to_csv(
            fname, index=False, mode='a', header=not os.path.exists(fname))
        print(f'saved {fname}')

    # plot profile
    df_snow_psd = pd.DataFrame.from_dict(psd_snow_dict)
    plot_psd_scatt_profile(
        df_snow_psd, args.path, band, hydro_label, ele=args.ele,
        x_var_list=psd_x_var_list, y_var_list=psd_y_var_list)


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
