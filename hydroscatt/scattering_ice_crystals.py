#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
scattering_ice_crystals
================================================

Scattering simulations of ice crystals


"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np
import pandas as pd

from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, ExponentialPSD
from pytmatrix import orientation, tmatrix_aux

from scattering_io import get_save_dir

from graph import plot_sp_scatt_quantities, plot_psd_scatt_quantities

from part_descrip import compute_ice_crystal_density, compute_axis_ratio_ic
from part_descrip import compute_ice_mass, compute_velocity_ic

from precip import compute_equi_rainfall_rate, compute_elwc, psd_ic_func

from scattering import compute_scattering_psd_tm, compute_scattering_sp_tm

from refractivity import refractive_index_ice, refractive_index_ice_particle
from refractivity import refractive_index_ic_func

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to ice crystals scattering simulations framework')

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
        default=-15.,
        help='air temperature. Default -15')

    parser.add_argument(
        '--ele', type=float,
        default=0.,
        help='elevation angle. Default 0 deg')

    args = parser.parse_args()

    print("====== ice crystals scattering simulation started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== ice crystals scattering simulation finished: ")

    compute_sp = True
    compute_psd = True
    with_subdirs = False
    hydro_type = 'ice_crystals'

    # parameters
    sp_var_list = [
        'sca_xsect_h', 'sca_xsect_v', 'refl_h', 'refl_v', 'ldr_h', 'ldr_v',
        'zdr', 'rho_hv', 'delta_hv', 'ext_xsect_h', 'ext_xsect_v', 'kdp',
        'A_h', 'A_v', 'Adp']

    sp_x_var_list = ['l']
    sp_y_var_list = [
        'sca_xsect', 'ext_xsect', 'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A', 'Adp']

    psd_var_list = [
        'refl_h', 'refl_v', 'ldr_h', 'ldr_v', 'zdr', 'rho_hv', 'delta_hv',
        'kdp', 'A_h', 'A_v', 'Adp']

    psd_x_var_list = ['refl_h', 'lwc', 'rr', 'D0']
    psd_y_var_list = [
        'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']

    ndgs = 2

    if args.temp < -25.:
        # Uppermost part of the cloud: small spherical ice particles dominant
        dominant_hydro = 'small_ice'
        length_min = 0.01  # length of the most elongated side (mm)
        length_max = 0.3
        step = 0.01

        area_ratio = 1.

        dens_coeff = 0.000916
        dens_exponent = 0.

        ar_coeff = 0.9
        ar_exponent = 1.
    elif -25. <= args.temp < -17.:
        # hexagonal plates dominant
        dominant_hydro = 'plates'
        length_min = 0.1
        length_max = 1.
        step = 0.1

        area_ratio = 0.83

        dens_coeff = 0.000916
        dens_exponent = 0.

        ar_coeff = 0.23
        ar_exponent = 0.778
    else:
        ndgs = 15

        # dendrites dominant
        dominant_hydro = 'dendrites'
        length_min = 0.1
        length_max = 2.4
        step = 0.1

        area_ratio = 0.28

        dens_coeff = 0.588e-3
        dens_exponent = -0.377

        ar_coeff = 0.0418
        ar_exponent = 0.377

    length = np.arange(length_min, length_max+step, step)
    num_points = length.size
    canting_angle = 15.

    # geometry = (theta0, theta, phi0, phi, alpha, beta)
    # geom_horiz_back = (90.0, 90.0, 0.0, 180.0, 0.0, 0.0) #horiz. backscatter
    # geom_horiz_forw = (90.0, 90.0, 0.0, 0.0, 0.0, 0.0) #horiz. forward scatter
    # geom_vert_back = (0.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. backscatter
    # geom_vert_forw = (180.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. forward scatter
    geom_back = (90.0-args.ele, 90.0+args.ele, 0.0, 180.0, 0.0, 0.0)
    geom_forw = (90.0-args.ele, 90.0-args.ele, 0.0, 0.0, 0.0, 0.0)

    print(hydro_type)
    print('dominant hydrometeor:', dominant_hydro)
    print('band:', args.band)
    print('temp:', args.temp)
    print('elevation angle', args.ele)

    if args.band == 'S':
        wavelength = tmatrix_aux.wl_S
    elif args.band == 'C':
        wavelength = tmatrix_aux.wl_C
    elif args.band == 'X':
        wavelength = tmatrix_aux.wl_X

    dens_ic = compute_ice_crystal_density(length, dens_coeff, dens_exponent)
    m_ice = refractive_index_ice(wavelength, args.temp)
    m_ic = refractive_index_ice_particle(m_ice, dens_ic)

    axis_ratio = compute_axis_ratio_ic(length, ar_coeff, ar_exponent)
    diam = length*np.power(axis_ratio, 1./3.)
    delta_d = diam-np.append(0., diam[:-1])
    mass = compute_ice_mass(diam, dens_ic)
    vel = compute_velocity_ic(length, mass, area_ratio)

    scatterer = Scatterer(wavelength=wavelength, ndgs=ndgs)
    scatterer.or_pdf = orientation.gaussian_pdf(canting_angle)
    scatterer.orient = orientation.orient_averaged_fixed
    scatterer.radius_type = Scatterer.RADIUS_MAXIMUM

    savedir = get_save_dir(
        args.path, hydro_type, args.band, args.temp, create_dir=True,
        with_subdirs=with_subdirs)

    if compute_sp:
        # single particle scattering
        single_part_dict = {'d': diam}
        single_part_dict.update({'l': length})
        for var in sp_var_list:
            single_part_dict.update({var: np.zeros(num_points)})

        for ind, l_part in enumerate(length):
            print(f'Computing point {ind} of {num_points} at D_max={l_part}')

            scatterer.m = m_ic[ind]
            scatterer.radius = l_part/2.
            scatterer.axis_ratio = 1./axis_ratio[ind]
            scatt_sp_dict = compute_scattering_sp_tm(
                scatterer, geom_back=geom_back, geom_forw=geom_forw,
                var_list=sp_var_list)

            for var in sp_var_list:
                single_part_dict[var][ind] = scatt_sp_dict[var]

        # save single particle results
        df_single_part = pd.DataFrame.from_dict(single_part_dict)
        fname = (
            f'{savedir}sp_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
            f'_ele{int(args.ele*100.):05d}_scattering.csv')
        df_single_part.to_csv(fname, index=False)
        print(f'saved {fname}')

        plot_sp_scatt_quantities(
            df_single_part, savedir, args.band, args.temp, hydro_type,
            ele=args.ele, x_var_list=sp_x_var_list, y_var_list=sp_y_var_list)

    if compute_psd:
        # Heymsfield et al (2013)
        # lambda = 1/mm
        lamb_lim = np.sort(
            np.array([0.34*np.exp(-0.083*args.temp),
                      1.53*np.exp(-0.053*args.temp)]))
        lam_vec_aux = np.arange(lamb_lim[0], lamb_lim[1], 0.1)
        n_d0 = lam_vec_aux.size

        # Nt in 1/m3
        nt_lim = np.sort(
            np.array([9.7e3*np.exp(-0.026*args.temp),
                      3.1e3*np.exp(-0.049*args.temp)]))
        nt_vec = np.arange(nt_lim[0], nt_lim[1], 10)
        n_nt = nt_vec.size

        n_dsd = n_d0*n_nt

        n0_vec = np.zeros(n_dsd)
        d0_vec = np.zeros(n_dsd)
        lamb_vec = np.zeros(n_dsd)
        for ind_lamb, lamb in enumerate(lam_vec_aux):
            for ind_nt, nt in enumerate(nt_vec):
                ind = ind_nt*n_d0+ind_lamb

                d0_vec[ind] = 3.67/lamb
                lamb_vec[ind] = lamb
                n0_vec[ind] = nt/np.sum(np.exp(-lamb*diam)*delta_d)

        print('nDSD:', n_dsd)
        print('lambda limits:', lamb_lim)
        print('Nt limits:', nt_lim)

        # DSD scattering properties
        scatterer.psd_integrator = PSDIntegrator()
        scatterer.psd_integrator.axis_ratio_func = (
            lambda length: 1.0/compute_axis_ratio_ic(
                    length, ar_coeff, ar_exponent))
        scatterer.psd_integrator.m_func = (
            lambda length: refractive_index_ic_func(
                length, dens_coeff, dens_exponent, m_ice))

        scatterer.psd_integrator.D_max = length[-1]
        scatterer.psd_integrator.num_points = num_points
        scatterer.psd_integrator.geometries = (geom_back, geom_forw)
        scatterer.psd_integrator.init_scatter_table(
            scatterer, angular_integration=True, verbose=True)

        psd_ic_dict = {
            'D0': d0_vec,
            'N0': n0_vec,
            'lwc': np.zeros(n_dsd),
            'rr': np.zeros(n_dsd),
        }
        for var in psd_var_list:
            psd_ic_dict.update({var: np.zeros(n_dsd)})

        for ind_dsd, (n0, lamb) in enumerate(zip(n0_vec, lamb_vec)):
            # lwc and rr computations are function of equivalent volume
            # diameter
            psd_d = ExponentialPSD(N0=n0, Lambda=lamb)
            psd_vals = psd_d(diam)
            psd_ic_dict['lwc'][ind_dsd] = compute_elwc(
                delta_d, mass, psd_vals)
            psd_ic_dict['rr'][ind_dsd] = compute_equi_rainfall_rate(
                delta_d, mass, psd_vals, vel)

            # scattering computation are function of maximum dimension
            scatterer.psd = lambda length: psd_ic_func(
                length, ar_coeff, ar_exponent, n0, lamb)
            scatt_psd_dict = compute_scattering_psd_tm(
                scatterer, geom_back=geom_back, geom_forw=geom_forw,
                var_list=psd_var_list)

            for var in psd_var_list:
                psd_ic_dict[var][ind_dsd] = scatt_psd_dict[var]

        # save DSD scattering results
        df_dsd = pd.DataFrame.from_dict(psd_ic_dict)
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
