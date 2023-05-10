#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
scattering_rain
================================================

Scattering simulations of rain

TODO:
- re-use computation of individual raindrops instead of repeating them

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
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, tmatrix_aux, refractive

from scattering_io import get_save_dir
from graph import plot_sp_scatt_quantities, plot_psd_scatt_quantities
from part_descrip import compute_velocity_rain
from precip import compute_lwc, compute_rainfall_rate
from scattering import compute_scattering_psd_tm, compute_scattering_sp_tm
from scattering import compute_scattering_canting_sp
from scattering import compute_scattering_canting_psd
from scattering import compute_angular_moments_analytical

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to rain scattering simulations framework')

    # keyword arguments
    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
        help='output data path')

    parser.add_argument(
        '--band', type=str,
        default='C',
        help='frequency band. Default C')

    parser.add_argument(
        '--temp', type=float,
        default=20.,
        help='air temperature. Default 20')

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

    print(f'====== rain scattering simulation started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== rain scattering simulation finished: ")

    if args.analytical_cant_angl == 1 and args.ele != 0:
        warn('The analytical computation of canting angle is only valid for '
             'elevation angles close to 0 deg')
        return

    with_subdirs = False
    hydro_type = 'rain'

    # parameters
    sp_var_list = [
        'sca_xsect_h', 'sca_xsect_v', 'refl_h', 'refl_v', 'ldr_h', 'ldr_v',
        'zdr', 'rho_hv', 'delta_hv', 'ext_xsect_h', 'ext_xsect_v', 'kdp',
        'A_h', 'A_v', 'Adp']

    sp_x_var_list = ['d']
    sp_y_var_list = [
        'sca_xsect', 'ext_xsect', 'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A', 'Adp']

    psd_var_list = [
        'refl_h', 'refl_v', 'ldr_h', 'ldr_v', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A_h', 'A_v', 'Adp']

    psd_x_var_list = ['refl_h', 'lwc', 'rr', 'D0']
    psd_y_var_list = [
        'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']

    diam_min = 0.1
    diam_max = 7.
    step = 0.1
    diam = np.arange(diam_min, diam_max+step, step)
    num_points = diam.size

    canting_angle = 10.

    # geometry = (theta0, theta, phi0, phi, alpha, beta)
    # geom_horiz_back = (90.0, 90.0, 0.0, 180.0, 0.0, 0.0) #horiz. backscatter
    # geom_horiz_forw = (90.0, 90.0, 0.0, 0.0, 0.0, 0.0) #horiz. forward scatter
    # geom_vert_back = (0.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. backscatter
    # geom_vert_forw = (180.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. forward scatter
    geom_back = (90.0-args.ele, 90.0+args.ele, 0.0, 180.0, 0.0, 0.0)
    geom_forw = (90.0-args.ele, 90.0-args.ele, 0.0, 0.0, 0.0, 0.0)

    print(hydro_type)
    print('band:', args.band)
    print('temp:', args.temp)
    print('elevation angle', args.ele)

    if args.band == 'S':
        wavelength = tmatrix_aux.wl_S
    elif args.band == 'C':
        wavelength = tmatrix_aux.wl_C
    elif args.band == 'X':
        wavelength = tmatrix_aux.wl_X

    if args.temp == 0.:
        m = refractive.m_w_0C[wavelength]
    elif args.temp == 10.:
        m = refractive.m_w_10C[wavelength]
    elif args.temp == 20.:
        m = refractive.m_w_20C[wavelength]

    scatterer = Scatterer(wavelength=wavelength, m=m)
    if args.analytical_cant_angl:
        scatterer.orient = orientation.orient_single
        ang_moments_dict = compute_angular_moments_analytical(canting_angle)
    else:
        scatterer.orient = orientation.orient_averaged_fixed
        scatterer.or_pdf = orientation.gaussian_pdf(canting_angle)

    savedir = get_save_dir(
        args.path, hydro_type, args.band, args.temp, create_dir=True,
        with_subdirs=with_subdirs)

    if args.compute_sp:
        # single particle scattering
        if args.analytical_cant_angl:
            fv180 = np.empty(num_points, dtype=complex)
            fh180 = np.empty(num_points, dtype=complex)
            fv0 = np.empty(num_points, dtype=complex)
            fh0 = np.empty(num_points, dtype=complex)
        else:
            single_part_dict = {'d': diam}
            for var in sp_var_list:
                single_part_dict.update({var: np.zeros(num_points)})

        for ind, d_part in enumerate(diam):
            print(f'Computing point {ind} at D={d_part}')

            scatterer.radius = d_part/2.
            scatterer.axis_ratio = 1.0/tmatrix_aux.dsr_thurai_2007(d_part)
            if args.analytical_cant_angl:
                scatterer.set_geometry(geom_back)
                s_mat = scatterer.get_S()
                fv180[ind] = 1e-3*s_mat[0][0]
                fh180[ind] = -1e-3*s_mat[1][1]
                scatterer.set_geometry(geom_forw)
                s_mat = scatterer.get_S()
                fv0[ind] = 1e-3*s_mat[0][0]
                fh0[ind] = 1e-3*s_mat[1][1]
            else:
                scatt_sp_dict = compute_scattering_sp_tm(
                    scatterer, geom_back=geom_back, geom_forw=geom_forw,
                    var_list=sp_var_list)

                for var in sp_var_list:
                    single_part_dict[var][ind] = scatt_sp_dict[var]

        if args.analytical_cant_angl:
            df_single_part = compute_scattering_canting_sp(
                wavelength, fv180, fh180, fv0, fh0, ang_moments_dict,
                var_list=sp_var_list)
            df_single_part['d'] = diam
        else:
            df_single_part = pd.DataFrame.from_dict(single_part_dict)

        # save single particle results
        fname = (
            f'{savedir}sp_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_single_part.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_sp_scatt_quantities(
            df_single_part, savedir, args.band, args.temp, hydro_type,
            ele=args.ele, x_var_list=sp_x_var_list, y_var_list=sp_y_var_list)

    if args.compute_psd:
        #    # Exponential Marshall-Palmer
        #    Nw_vec = np.array([8e3])
        #    nNw = Nw_vec.size
        #
        #    mu_vec = np.array([0])
        #    nmu = mu_vec.size
        #
        #    rmax = 300  # mm/h
        #    rmin = 0.5
        #    step = 0.5
        #
        #    D0_vec = 3.67/(4.1*(np.arange(rmin, rmax+step, step))**(-0.21))
        #    nD0 = D0_vec.size

        # Generic gamma distribution
        rmin = 1
        rmax = 4
        step = 0.1
        nw_vec_aux = np.power(10, np.arange(rmin, rmax+step, step))
        n_nw = nw_vec_aux.size

        rmin = -1
        rmax = 8
        step = 1
        mu_vec_aux = np.arange(rmin, rmax+step, step)
        n_mu = mu_vec_aux.size

        rmax = 3
        rmin = 0.1
        step = 0.1
        d0_vec_aux = np.arange(rmin, rmax+step, step)
        n_d0 = d0_vec_aux.size

        n_dsd = n_d0*n_mu*n_nw

        d0_vec = np.zeros(n_dsd)
        nw_vec = np.zeros(n_dsd)
        mu_vec = np.zeros(n_dsd)
        for ind_d0, d0 in enumerate(d0_vec_aux):
            for ind_nw, nw in enumerate(nw_vec_aux):
                for ind_mu, mu in enumerate(mu_vec_aux):
                    ind = ind_mu*n_d0*n_nw+ind_nw*n_d0+ind_d0

                    d0_vec[ind] = d0
                    nw_vec[ind] = nw
                    mu_vec[ind] = mu

        print('nDSD:', n_dsd)

        delta_d = diam-np.append(0., diam[:-1])
        vel = compute_velocity_rain(diam, rho=1.22)

        if not args.analytical_cant_angl:
            # DSD scattering properties
            scatterer.psd_integrator = PSDIntegrator()
            scatterer.psd_integrator.axis_ratio_func = (
                lambda diam: 1.0/tmatrix_aux.dsr_thurai_2007(diam))
            scatterer.psd_integrator.D_max = diam_max
            scatterer.psd_integrator.num_points = num_points
            scatterer.psd_integrator.geometries = (geom_back, geom_forw)
            scatterer.psd_integrator.init_scatter_table(
                scatterer, angular_integration=True, verbose=True)

        dsd_rain_dict = {
            'D0': d0_vec,
            'Nw': nw_vec,
            'mu': mu_vec,
            'lwc': np.zeros(n_dsd),
            'rr': np.zeros(n_dsd),
        }
        for var in psd_var_list:
            dsd_rain_dict.update({var: np.zeros(n_dsd)})

        for ind_dsd in range(n_dsd):
            psd = GammaPSD(
                D0=d0_vec[ind_dsd], Nw=nw_vec[ind_dsd], mu=mu_vec[ind_dsd],
                D_max=diam_max)
            psd_vals = psd(diam)
            dsd_rain_dict['lwc'][ind_dsd] = compute_lwc(
                diam, delta_d, psd_vals)
            dsd_rain_dict['rr'][ind_dsd] = compute_rainfall_rate(
                diam, delta_d, psd_vals, vel)

            if args.analytical_cant_angl:
                scatt_psd_dict = compute_scattering_canting_psd(
                    wavelength, fv180, fh180, fv0, fh0, ang_moments_dict,
                    delta_d, psd_vals, var_list=psd_var_list)
            else:
                scatterer.psd = psd
                scatt_psd_dict = compute_scattering_psd_tm(
                    scatterer, geom_back=geom_back, geom_forw=geom_forw,
                    var_list=psd_var_list)

            for var in psd_var_list:
                dsd_rain_dict[var][ind_dsd] = scatt_psd_dict[var]

        # save DSD scattering results
        df_dsd = pd.DataFrame.from_dict(dsd_rain_dict)
        fname = (
            f'{savedir}psd_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_dsd.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_psd_scatt_quantities(
            df_dsd, savedir, args.band, args.temp, hydro_type, ele=args.ele,
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
