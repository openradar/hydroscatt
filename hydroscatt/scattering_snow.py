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
from warnings import warn

import numpy as np
import pandas as pd

from pytmatrix.psd import ExponentialPSD

from scattering_io import get_save_dir, read_melting_hydro_part_model
from scattering_io import read_scatt_double_layer, read_wavelength_file

from precip import compute_elwc, compute_equi_rainfall_rate

from refractivity import wavelength_to_band

from scattering import compute_angular_moments, compute_scattering_canting_sp
from scattering import compute_scattering_canting_psd
from scattering import compute_angular_moments_analytical

from graph import plot_sp_scatt_quantities, plot_psd_scatt_quantities

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
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='input data path')

    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='data path')

    parser.add_argument(
        '--model_file', type=str,
        default='sp_melting_snow_C_0000_model_part.csv',
        help='particle model file')

    parser.add_argument(
        '--scatt_file', type=str,
        default='sp_melting_snow_C_0000_tmat.out',
        help='scattering matrix elements file')

    parser.add_argument(
        '--freq_file', type=str,
        default='sp_melting_snow_C_0000_freq.inp',
        help='frequency file')

    parser.add_argument(
        '--ele', type=float,
        default=0.,
        help='elevation angle. Default 0 deg')

    parser.add_argument(
        '--compute_sp', type=int,
        default=1,
        help='whether to compute single particle scattering. Default 1')

    parser.add_argument(
        '--compute_psd', type=int,
        default=1,
        help='whether to compute PSD scattering. Default 1')

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

    with_subdirs = False
    hydro_type = 'melting_snow'

    # parameters
    sp_var_list = [
        'sca_xsect_h', 'sca_xsect_v', 'refl_h', 'refl_v', 'ldr_h', 'ldr_v',
        'zdr', 'rho_hv', 'delta_hv', 'ext_xsect_h', 'ext_xsect_v', 'kdp',
        'A_h', 'A_v', 'Adp']

    sp_x_var_list = ['d_ext']
    sp_y_var_list = [
        'sca_xsect', 'ext_xsect', 'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A', 'Adp']

    psd_var_list = [
        'refl_h', 'refl_v', 'ldr_h', 'ldr_v', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A_h', 'A_v', 'Adp']

    psd_x_var_list = ['refl_h', 'lwc', 'rr', 'D0']
    psd_y_var_list = [
        'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']

    df_model, temp = read_melting_hydro_part_model(
        f'{args.input_path}{args.model_file}')
    df_scatt = read_scatt_double_layer(f'{args.input_path}{args.scatt_file}')
    wavelength = read_wavelength_file(f'{args.input_path}{args.freq_file}')

    band = wavelength_to_band(wavelength)

    print(hydro_type)
    print('temp', temp)
    print('band', band)
    print('elevation angle', args.ele)

    print('computing angular moments ...')
    canting_angle_snow = 40.+df_model['fmw'].values*(10-40)
    if args.analytical_cant_angl:
        ang_moments_dict = compute_angular_moments_analytical(
            canting_angle_snow)
    else:
        ang_moments_dict = compute_angular_moments(
            canting_angle_snow, ele=args.ele)

    savedir = get_save_dir(
        args.path, hydro_type, band, temp, create_dir=True,
        with_subdirs=with_subdirs)

    if args.compute_sp:
        print(
            'computing single scattering parameters of melting snowflakes ...')
        df_snow_sp = compute_scattering_canting_sp(
            wavelength, df_scatt['fv180'].values, df_scatt['fh180'].values,
            df_scatt['fv0'].values, df_scatt['fh0'].values, ang_moments_dict,
            var_list=sp_var_list)

        df_snow_sp['d_init'] = df_model['d_init']
        df_snow_sp['d_ext'] = df_model['d_ext']
        fname = (
            f'{savedir}sp_{hydro_type}_{band}_{int(temp*100.):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_snow_sp.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_sp_scatt_quantities(
            df_snow_sp, savedir, band, temp, hydro_type, ele=args.ele,
            x_var_list=sp_x_var_list, y_var_list=sp_y_var_list)

    if args.compute_psd:
        print('computing PSD scattering parameters of melting snowflakes ...')

        # Exponential distribution (Brandes et al (2007))
        # D0 = 3.67/lambda
        rmin = 1
        rmax = 6.
        step = 0.1
        lamb_vec_aux = np.arange(rmin, rmax+step, step)
        n_lamb = lamb_vec_aux.size

        rmin = 3
        rmax = 4.5
        step = 0.01
        nw_vec_aux = np.power(10., np.arange(rmin, rmax+step, step))
        n_nw = nw_vec_aux.size

        n_dsd = n_lamb*n_nw

        nw_vec = np.zeros(n_dsd)
        d0_vec = np.zeros(n_dsd)
        lamb_vec = np.zeros(n_dsd)
        for ilamb, lamb in enumerate(lamb_vec_aux):
            for inw, nw in enumerate(nw_vec_aux):
                ind = ilamb+inw*n_lamb

                nw_vec[ind] = nw
                lamb_vec[ind] = lamb
                d0_vec[ind] = 3.67/lamb

        print('nDSD', n_dsd)

        delta_d_snow = (
            df_model['d_ext'].values
            - np.append(0, df_model['d_ext'].values[:-1]))
        vel_snow = df_model['vel'].values

        psd_snow_dict = {
            'D0': d0_vec,
            'Nw': nw_vec,
            'lwc': np.zeros(n_dsd),
            'rr': np.zeros(n_dsd)
        }
        for var in psd_var_list:
            psd_snow_dict.update({var: np.zeros(n_dsd)})

        for ind_dsd in range(n_dsd):
            psd_snow = ExponentialPSD(
                N0=nw_vec[ind_dsd], Lambda=lamb_vec[ind_dsd],
                D_max=df_model['d_init'].values[-1])  # should be defined by D0
            psd_vals_snow = psd_snow(df_model['d_init'].values)

            psd_snow_dict['lwc'][ind_dsd] = compute_elwc(
                delta_d_snow, df_model['mass'].values, psd_vals_snow)
            psd_snow_dict['rr'][ind_dsd] = compute_equi_rainfall_rate(
                delta_d_snow, df_model['mass'].values, psd_vals_snow,
                vel_snow)

            snow_psd_dict = compute_scattering_canting_psd(
                wavelength, df_scatt['fv180'].values, df_scatt['fh180'].values,
                df_scatt['fv0'].values, df_scatt['fh0'].values,
                ang_moments_dict, delta_d_snow, psd_vals_snow,
                var_list=psd_var_list)

            for var in psd_var_list:
                psd_snow_dict[var][ind_dsd] = snow_psd_dict[var]

        # save DSD scattering results
        df_snow_psd = pd.DataFrame.from_dict(psd_snow_dict)
        fname = (
            f'{savedir}psd_{hydro_type}_{band}_{int(temp*100.):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_snow_psd.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_psd_scatt_quantities(
            df_snow_psd, savedir, band, temp, hydro_type, ele=args.ele,
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
