#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
scattering_hail
================================================

modelling of scattering properties of hail

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import glob
import argparse
import atexit
import os
from warnings import warn
from copy import deepcopy

import numpy as np
import pandas as pd

from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator
from pytmatrix import orientation, tmatrix_aux

from scattering_io import get_save_dir, read_melting_hydro_part_model
from scattering_io import read_scatt_double_layer, read_wavelength_file

from part_descrip import compute_velocity_rain

from precip import compute_lwc, compute_elwc, compute_rainfall_rate
from precip import compute_equi_rainfall_rate, compute_dsd_shed_water
from precip import psd_hail, compute_dsd_breakup_water

from refractivity import wavelength_to_band, refractive_index_water

from scattering import compute_angular_moments, compute_scattering_canting_psd
from scattering import compute_scattering_psd_tm, compute_scattering_mixture
from scattering import compute_angular_moments_analytical

from graph import plot_psd_scatt_profile

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to hail scattering simulations framework')

    # keyword arguments
    parser.add_argument(
        '--input_path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/profile/',
        help='input data path')

    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/profile/',
        help='data path')

    parser.add_argument(
        '--band', type=str,
        default='C',
        help='frequency band. Default C')

    parser.add_argument(
        '--hail_size', type=str,
        default='large',
        help='hail size. Can be large, moderate or small. Default large')

    parser.add_argument(
        '--ele', type=float,
        default=0.,
        help='elevation angle. Default 0 deg')

    parser.add_argument(
        '--analytical_cant_angl', type=int,
        default=1,
        help='If 1 the canting angle will be computed analytically. Default 1')

    args = parser.parse_args()

    print(
        f"====== hail scattering simulation started: "
        f"{datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')}")
    atexit.register(_print_end_msg,
                    "====== hail scattering simulation finished: ")

    if args.analytical_cant_angl == 1 and args.ele != 0:
        warn('The analytical computation of canting angle is only valid for '
             'elevation angles close to 0 deg')
        return

    hydro_type = 'melting_hail'

    iso0_height = 4000.
    lapse_rate = 6.5

    if args.hail_size == 'large':
        # PSD for large hail
        lamb_h = 0.27
        coeff_nh = 800.
        d_max_h = 35.
        hydro_label = 'large_melting_hail'
    elif args.hail_size == 'moderate':
        # PSD for moderate hail
        lamb_h = 0.42
        coeff_nh = 400.
        d_max_h = 24.
        hydro_label = 'moderate_melting_hail'
    else:
        # PSD for small hail
        lamb_h = 0.99
        coeff_nh = 200.
        d_max_h = 14.
        hydro_label = 'small_melting_hail'

    # parameters for rain
    canting_angle_rain = 10
    d_drop_min = 0.1
    d_shed_max = 4.5
    step = 0.1
    d_shed = np.arange(d_drop_min, d_shed_max+step, step)

    # parameters
    psd_var_list = [
        'refl_h', 'refl_v', 'ldr_h', 'ldr_v', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A_h', 'A_v', 'Adp']

    psd_x_var_list = [
        'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']
    psd_y_var_list = ['temp']

    # geometry = (theta0, theta, phi0, phi, alpha, beta)
    # geom_horiz_back = (90.0, 90.0, 0.0, 180.0, 0.0, 0.0) #horiz. backscatter
    # geom_horiz_forw = (90.0, 90.0, 0.0, 0.0, 0.0, 0.0) #horiz. forward scatter
    # geom_vert_back = (0.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. backscatter
    # geom_vert_forw = (180.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. forward scatter
    geom_back = (90.0-args.ele, 90.0+args.ele, 0.0, 180.0, 0.0, 0.0)
    geom_forw = (90.0-args.ele, 90.0-args.ele, 0.0, 0.0, 0.0, 0.0)

    print(hydro_type)
    print(hydro_label)
    print('band', args.band)
    print('elevation angle', args.ele)

    flist_model = glob.glob(
        f'{args.input_path}sp_*{hydro_type}_*_model_part.csv')
    flist_scatt = glob.glob(
        f'{args.input_path}sp_*{hydro_type}_{args.band}_*_tmat.out')
    if not flist_model or not flist_scatt:
        if not flist_model:
            warn(f'no model file at '
                 f'{args.input_path}sp_*{hydro_type}_*_model_part.csv')
        if not flist_scatt:
            warn(f'no scattering file at '
                 f'{args.input_path}sp_*{hydro_type}_{args.band}'
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
    psd_hail_dict = {
        'temp': np.zeros(ntemp),
        'alt': np.zeros(ntemp),
        'lwc': np.zeros(ntemp),
        'rr': np.zeros(ntemp)
    }
    for var in psd_var_list:
        psd_hail_dict.update({var: np.zeros(ntemp)})

    if args.analytical_cant_angl:
        ang_moments_rain_dict = compute_angular_moments_analytical(
            canting_angle_rain)

    # output file name
    savedir = get_save_dir(
        args.path, hydro_type, args.band, None, create_dir=True,
        with_subdirs=False)

    fname = (
        f'{savedir}psd_profile_{hydro_label}_{band}'
        f'_ele{int(args.ele*100.):05d}_scattering.csv')
    for ind_temp, (model_file, scatt_file) in enumerate(zip(
            flist_model, flist_scatt)):
        df_model, temp = read_melting_hydro_part_model(
            model_file, d_max=d_max_h)
        n_hs = df_model.shape[0]
        df_scatt = read_scatt_double_layer(scatt_file, nrows=n_hs)

        alt = iso0_height - temp/lapse_rate*1000.

        psd_hail_dict['temp'][ind_temp] = temp
        psd_hail_dict['alt'][ind_temp] = alt
        print('temp', temp, 'alt', alt)
        print('computing angular moments ...')
        canting_angle_hail = 40.+df_model['fmw'].values*(10-40)
        if args.analytical_cant_angl:
            ang_moments_dict = compute_angular_moments_analytical(
                canting_angle_hail)
        else:
            ang_moments_dict = compute_angular_moments(
                canting_angle_hail, ele=args.ele)

        print('computing PSD scattering parameters of melting hailstones ...')

        shed_water_present = False
        if np.sum(df_model['shed_water_mass'].values) > 0.:
            shed_water_present = True

        breakup_water_present = False
        d_breakup_max = 0.
        if np.sum(df_model['prob_break'].values) > 0.:
            d_breakup = df_model['d_ext'][np.isclose(df_model['fmw'], 1.)]
            d_breakup_max = int(10.*np.max(d_breakup))/10.
            breakup_water_present = True

        # single raindrop scattering
        if shed_water_present or breakup_water_present:
            m_water = refractive_index_water(wavelength, temp, sal=0.)
            d_drop_max = np.max((d_shed_max, d_breakup_max))
            d_drop = np.arange(d_drop_min, d_drop_max+step, step)
            delta_d_drop = d_drop-np.append(0, d_drop[:-1])

            scatterer = Scatterer(wavelength=wavelength, m=m_water)
            if args.analytical_cant_angl:
                fv180 = np.empty(d_drop.size, dtype=complex)
                fh180 = np.empty(d_drop.size, dtype=complex)
                fv0 = np.empty(d_drop.size, dtype=complex)
                fh0 = np.empty(d_drop.size, dtype=complex)
                scatterer.orient = orientation.orient_single
                for ind, d_part in enumerate(d_drop):
                    print(f'Computing point {ind} at D={d_part}')
                    scatterer.radius = d_part/2.
                    scatterer.axis_ratio = 1.0/tmatrix_aux.dsr_thurai_2007(
                        d_part)
                    scatterer.set_geometry(geom_back)
                    s_mat = scatterer.get_S()
                    fv180[ind] = 1e-3*s_mat[0][0]
                    fh180[ind] = -1e-3*s_mat[1][1]
                    scatterer.set_geometry(geom_forw)
                    s_mat = scatterer.get_S()
                    fv0[ind] = 1e-3*s_mat[0][0]
                    fh0[ind] = 1e-3*s_mat[1][1]
            else:
                scatterer.or_pdf = orientation.gaussian_pdf(canting_angle_rain)
                scatterer.orient = orientation.orient_averaged_fixed
                scatterer.psd_integrator = PSDIntegrator()
                scatterer.psd_integrator.axis_ratio_func = (
                    lambda d_drop: 1.0/tmatrix_aux.dsr_thurai_2007(d_drop))
                scatterer.psd_integrator.D_max = d_drop_max
                scatterer.psd_integrator.num_points = d_drop.size
                scatterer.psd_integrator.geometries = (geom_back, geom_forw)
                scatterer.psd_integrator.init_scatter_table(
                    scatterer, angular_integration=True, verbose=True)

        vel_hail = df_model['vel'].values

        if np.isclose(alt, iso0_height):
            vel_hail0 = deepcopy(vel_hail)

        psd_hail_no_break = psd_hail(
            df_model['d_init'].values, lamb_h=lamb_h, coeff_nh=coeff_nh,
            d_max_h=df_model['d_init'].values[-1])

        # conservation of the flux
        psd_hail_no_break *= (vel_hail0/vel_hail)

        psd_break_vals = psd_hail_no_break*df_model['prob_break'].values
        psd_vals_hail = psd_hail_no_break-psd_break_vals

        delta_d_hail = (
            df_model['d_init'].values
            - np.append(0, df_model['d_init'].values[:-1]))

        psd_hail_dict['lwc'][ind_temp] = compute_elwc(
            delta_d_hail, df_model['hail_mass'].values, psd_vals_hail)
        psd_hail_dict['rr'][ind_temp] = compute_equi_rainfall_rate(
            delta_d_hail, df_model['hail_mass'].values, psd_vals_hail,
            vel_hail)

        hail_psd_dict = compute_scattering_canting_psd(
            wavelength, df_scatt['fv180'].values, df_scatt['fh180'].values,
            df_scatt['fv0'].values, df_scatt['fh0'].values,
            ang_moments_dict, delta_d_hail, psd_vals_hail,
            var_list=psd_var_list)

        if shed_water_present or breakup_water_present:
            psd_rain_vals = np.zeros(d_drop.size)
            if shed_water_present:
                psd_shed_func = compute_dsd_shed_water(
                    df_model['d_init'].values, psd_hail_no_break, d_shed,
                    df_model['shed_water_mass'].values)
                psd_rain_vals += psd_shed_func(d_drop)

                if not args.analytical_cant_angl:
                    scatterer.psd = psd_shed_func
                    rain_psd_dict = compute_scattering_psd_tm(
                        scatterer, geom_back=geom_back, geom_forw=geom_forw,
                        var_list=psd_var_list)

                    hail_psd_dict = compute_scattering_mixture(
                        hail_psd_dict, rain_psd_dict, var_list=psd_var_list)

            if breakup_water_present:
                psd_break_func = compute_dsd_breakup_water(
                    df_model['d_init'].values, psd_break_vals, d_breakup,
                    df_model['hail_mass'].values*df_model['fmw'].values)
                psd_rain_vals += psd_break_func(d_drop)

                if not args.analytical_cant_angl:
                    scatterer.psd = psd_break_func
                    rain_psd_dict = compute_scattering_psd_tm(
                        scatterer, geom_back=geom_back, geom_forw=geom_forw,
                        var_list=psd_var_list)

                    hail_psd_dict = compute_scattering_mixture(
                        hail_psd_dict, rain_psd_dict, var_list=psd_var_list)

            if args.analytical_cant_angl:
                rain_psd_dict = compute_scattering_canting_psd(
                    wavelength, fv180, fh180, fv0, fh0, ang_moments_rain_dict,
                    delta_d_drop, psd_rain_vals, var_list=psd_var_list)

                hail_psd_dict = compute_scattering_mixture(
                    hail_psd_dict, rain_psd_dict, var_list=psd_var_list)

            vel_rain = compute_velocity_rain(d_drop)
            lwc_rain = compute_lwc(d_drop, delta_d_drop, psd_rain_vals)
            rr_rain = compute_rainfall_rate(
                d_drop, delta_d_drop, psd_rain_vals, vel_rain)

            psd_hail_dict['lwc'][ind_temp] += lwc_rain
            psd_hail_dict['rr'][ind_temp] += rr_rain

        for var in psd_var_list:
            psd_hail_dict[var][ind_temp] = hail_psd_dict[var]

        # save DSD scattering results
        psd_hail_dict_aux = {}
        for key, value in psd_hail_dict.items():
            psd_hail_dict_aux.update({key: [value[ind_temp]]})
        df_hail_psd_aux = pd.DataFrame.from_dict(psd_hail_dict_aux)
        df_hail_psd_aux.to_csv(
            fname, index=False, mode='a', header=not os.path.exists(fname))
        print(f'saved {fname}')

    # plot profile
    df_hail_psd = pd.DataFrame.from_dict(psd_hail_dict)
    plot_psd_scatt_profile(
        df_hail_psd, savedir, band, hydro_label, ele=args.ele,
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
