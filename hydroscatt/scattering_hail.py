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
import argparse
import atexit
from warnings import warn

import numpy as np
import pandas as pd

from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, ExponentialPSD
from pytmatrix import orientation, tmatrix_aux

from scattering_io import get_save_dir, read_melting_hydro_part_model
from scattering_io import read_scatt_double_layer, read_wavelength_file

from part_descrip import compute_velocity_rain

from precip import compute_lwc, compute_elwc, compute_rainfall_rate
from precip import compute_equi_rainfall_rate, compute_dsd_shed_water
from precip import compute_dsd_breakup_water

from refractivity import wavelength_to_band, refractive_index_water

from scattering import compute_angular_moments, compute_scattering_canting_sp
from scattering import compute_scattering_canting_psd
from scattering import compute_scattering_psd_tm, compute_scattering_mixture
from scattering import compute_angular_moments_analytical

from graph import plot_sp_scatt_quantities, plot_psd_scatt_quantities

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
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='input data path')

    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='data path')

    parser.add_argument(
        '--model_file', type=str,
        default='sp_melting_hail_0000_model_part.csv',
        help='particle model file')

    parser.add_argument(
        '--scatt_file', type=str,
        default='sp_melting_hail_C_0000_tmat.out',
        help='scattering matrix elements file')

    parser.add_argument(
        '--freq_file', type=str,
        default='sp_melting_hail_C_0000_freq.inp',
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

    print(f'====== hail scattering simulation started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== hail scattering simulation finished: ")

    if args.analytical_cant_angl == 1 and args.ele != 0:
        warn('The analytical computation of canting angle is only valid for '
             'elevation angles close to 0 deg')
        return

    with_subdirs = False
    hydro_type = 'melting_hail'

    canting_angle_rain = 10
    d_drop_min = 0.1
    d_shed_max = 4.5
    d_drop_step = 0.1

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

    # geometry = (theta0, theta, phi0, phi, alpha, beta)
    # geom_horiz_back = (90.0, 90.0, 0.0, 180.0, 0.0, 0.0) #horiz. backscatter
    # geom_horiz_forw = (90.0, 90.0, 0.0, 0.0, 0.0, 0.0) #horiz. forward scatter
    # geom_vert_back = (0.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. backscatter
    # geom_vert_forw = (180.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. forward scatter
    geom_back = (90.0-args.ele, 90.0+args.ele, 0.0, 180.0, 0.0, 0.0)
    geom_forw = (90.0-args.ele, 90.0-args.ele, 0.0, 0.0, 0.0, 0.0)

    print(hydro_type)
    print('temp', temp)
    print('band', band)
    print('elevation angle', args.ele)

    print('computing angular moments ...')
    canting_angle_hail = 40.+df_model['fmw'].values*(10-40)
    if args.analytical_cant_angl:
        ang_moments_dict = compute_angular_moments_analytical(
            canting_angle_hail)
        ang_moments_rain_dict = compute_angular_moments_analytical(
            canting_angle_rain)
    else:
        ang_moments_dict = compute_angular_moments(
            canting_angle_hail, ele=args.ele)

    savedir = get_save_dir(
        args.path, hydro_type, band, temp, create_dir=True,
        with_subdirs=with_subdirs)

    if args.compute_sp:
        print(
            'computing single scattering parameters of melting hailstones ...')
        df_hail_sp = compute_scattering_canting_sp(
            wavelength, df_scatt['fv180'].values, df_scatt['fh180'].values,
            df_scatt['fv0'].values, df_scatt['fh0'].values, ang_moments_dict,
            var_list=sp_var_list)

        df_hail_sp['d_init'] = df_model['d_init']
        df_hail_sp['d_ext'] = df_model['d_ext']

        fname = (
            f'{savedir}sp_{hydro_type}_{band}_{int(temp*100.):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_hail_sp.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_sp_scatt_quantities(
            df_hail_sp, savedir, band, temp, hydro_type, ele=args.ele,
            x_var_list=sp_x_var_list, y_var_list=sp_y_var_list)

    if args.compute_psd:
        print('computing PSD scattering parameters of melting hailstones ...')

        # exponential distribution (Cheng et al, 1985)
        # N(D) = N0*exp(-lambda*D)
        # N0 = C*lambda^4.11
        # D0 = 3.67/lambda
        rmin = 0.1
        rmax = 1.
        step = 0.01
        lamb_vec_aux = np.arange(rmin, rmax+step, step)
        n_lamb = lamb_vec_aux.size

        rmin = 60
        rmax = 300
        step = 2
        c_par_vec = np.arange(rmin, rmax+step, step)
        n_c_par = c_par_vec.size

        n_dsd = n_lamb*n_c_par

        n0_vec = np.zeros(n_dsd)
        d0_vec = np.zeros(n_dsd)
        lamb_vec = np.zeros(n_dsd)
        for ilamb, lamb in enumerate(lamb_vec_aux):
            for ic_par, c_par in enumerate(c_par_vec):
                ind = ilamb+ic_par*n_lamb

                n0_vec[ind] = c_par*np.power(lamb, 4.11)
                lamb_vec[ind] = lamb
                d0_vec[ind] = 3.67/lamb

        print('nDSD', n_dsd)

        shed_water_present = False
        if np.sum(df_model['shed_water_mass'].values) > 0:
            d_shed = np.arange(d_drop_min, d_shed_max+d_drop_step, d_drop_step)
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
            d_drop = np.arange(d_drop_min, d_drop_max+d_drop_step, d_drop_step)
            delta_d_drop = d_drop-np.append(0, d_drop[:-1])
            vel_rain = compute_velocity_rain(d_drop)

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

        delta_d_hail = (
            df_model['d_init'].values
            - np.append(0, df_model['d_init'].values[:-1]))
        vel_hail = df_model['vel'].values

        psd_hail_dict = {
            'D0': d0_vec,
            'N0': n0_vec,
            'lwc': np.zeros(n_dsd),
            'rr': np.zeros(n_dsd)
        }
        for var in psd_var_list:
            psd_hail_dict.update({var: np.zeros(n_dsd)})

        for ind_dsd in range(n_dsd):
            psd_hail_func = ExponentialPSD(
                N0=n0_vec[ind_dsd], Lambda=lamb_vec[ind_dsd],
                D_max=df_model['d_init'].values[-1])
            psd_hail_no_break = psd_hail_func(df_model['d_init'].values)
            # conservation of the flux
            # psd_hail_no_break *= (vel_hail0/vel_hail)

            psd_break_vals = psd_hail_no_break*df_model['prob_break'].values
            psd_vals_hail = psd_hail_no_break-psd_break_vals

            psd_hail_dict['lwc'][ind_dsd] = compute_elwc(
                delta_d_hail, df_model['hail_mass'].values, psd_vals_hail)
            psd_hail_dict['rr'][ind_dsd] = compute_equi_rainfall_rate(
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
                            scatterer, geom_back=geom_back,
                            geom_forw=geom_forw, var_list=psd_var_list)

                        hail_psd_dict = compute_scattering_mixture(
                            hail_psd_dict, rain_psd_dict,
                            var_list=psd_var_list)

                if breakup_water_present:
                    psd_break_func = compute_dsd_breakup_water(
                        df_model['d_init'].values, psd_break_vals, d_breakup,
                        df_model['hail_mass'].values*df_model['fmw'].values)
                    psd_rain_vals += psd_break_func(d_drop)

                    if not args.analytical_cant_angl:
                        scatterer.psd = psd_break_func
                        rain_psd_dict = compute_scattering_psd_tm(
                            scatterer, geom_back=geom_back,
                            geom_forw=geom_forw, var_list=psd_var_list)

                        hail_psd_dict = compute_scattering_mixture(
                            hail_psd_dict, rain_psd_dict,
                            var_list=psd_var_list)

                if args.analytical_cant_angl:
                    rain_psd_dict = compute_scattering_canting_psd(
                        wavelength, fv180, fh180, fv0, fh0,
                        ang_moments_rain_dict, delta_d_drop, psd_rain_vals,
                        var_list=psd_var_list)

                    hail_psd_dict = compute_scattering_mixture(
                        hail_psd_dict, rain_psd_dict, var_list=psd_var_list)

                lwc_rain = compute_lwc(d_drop, delta_d_drop, psd_rain_vals)
                rr_rain = compute_rainfall_rate(
                    d_drop, delta_d_drop, psd_rain_vals, vel_rain)

                psd_hail_dict['lwc'][ind_dsd] += lwc_rain
                psd_hail_dict['rr'][ind_dsd] += rr_rain

            for var in psd_var_list:
                psd_hail_dict[var][ind_dsd] = hail_psd_dict[var]

        # save DSD scattering results
        df_hail_psd = pd.DataFrame.from_dict(psd_hail_dict)
        fname = (
            f'{savedir}psd_{hydro_type}_{band}_{int(temp*100.):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_hail_psd.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_psd_scatt_quantities(
            df_hail_psd, savedir, band, temp, hydro_type, ele=args.ele,
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
