#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
model_hail_profile
================================================

modelling of the characteristics of melting hail

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
from copy import deepcopy

import numpy as np

from pytmatrix import tmatrix_aux

from part_descrip import compute_ice_mass, compute_hailstone_density
from part_descrip import compute_hs_air_vol, compute_hail_melting
from part_descrip import hailstone_change, compute_axis_ratio_dh
from part_descrip import compute_axis_ratio_mh

from atmos import compute_air_pressure, compute_air_density
from atmos import compute_dynamic_viscosity

from precip import compute_breakup_prob

from refractivity import refractive_index_water
from refractivity import refractive_index_ice
from refractivity import refractive_index_ice_particle
from refractivity import refractive_index_melting_hail_core

from scattering_io import get_save_dir, write_melting_hydro_scatt_model
from scattering_io import write_melting_hail_part_model, write_wavelength_file

from scattering import compute_tmatrix_2l

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
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/profile/',
        help='output data path')

    parser.add_argument(
        '--tm2l_dir', type=str,
        default='/home/mdso/figuerasiventuraj/hydroscatt/hydroscatt/f_2l_tmatrix/',
        help='Fortran 2-layer T-matrix directory')

    parser.add_argument(
        '--band', type=str,
        default='C',
        help='frequency band. Default C')

    parser.add_argument(
        '--temp', type=float,
        default=0.,
        help='air temperature. Default 0')

    parser.add_argument(
        '--compute_scatt', type=int,
        default=1,
        help='Whether to compute the scattering matrix. Default 1 (True)')

    parser.add_argument(
        '--variable_density', type=int,
        default=1,
        help=(
            'Whether to use variable density hail or fixed high density hail'
            'Default 1 (True)'))

    args = parser.parse_args()

    print(f'====== hail modelling started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== hail modelling finished: ")

    with_subdirs = False

    hydro_type = 'melting_hail'

    iso0_height = 4000.

    temp_min = 0.
    lapse_rate = 6.5  # deg C/km
    delta_h = 10.  # m
    step = delta_h*lapse_rate/1000.

    temp_vec = np.arange(temp_min, args.temp+step, step)
    alt_vec = iso0_height - temp_vec/lapse_rate*1000.
    p_air_vec = compute_air_pressure(alt_vec)
    dens_air_vec = compute_air_density(alt_vec)
    visc_vec = compute_dynamic_viscosity(temp_vec)

    d_max = 35.
    d_min = 0.1
    step = 0.1
    diam = np.arange(d_min, d_max+step, step)

    if args.band == 'S':
        wavelength = tmatrix_aux.wl_S
    elif args.band == 'C':
        wavelength = tmatrix_aux.wl_C
    elif args.band == 'X':
        wavelength = tmatrix_aux.wl_X

    print(hydro_type)
    print('band:', args.band)
    print('temp:', args.temp)

    dens_hail_init = compute_hailstone_density(
        diam, variable_density=args.variable_density)
    mass_hail_init = compute_ice_mass(diam, dens_ice=dens_hail_init)
    vol_soakable = compute_hs_air_vol(diam, mass_hail_init, dens_hail_init)

    d_ext = deepcopy(diam)
    d_int = deepcopy(diam)
    fmw = np.zeros(diam.size)
    alpha = np.zeros(diam.size)+np.nan
    mass_hail = deepcopy(mass_hail_init)

    savedir = get_save_dir(
        args.path, hydro_type, args.band, args.temp, create_dir=True,
        with_subdirs=with_subdirs)

    fname_freq = f'{savedir}sp_{hydro_type}_{args.band}_freq.inp'
    fname_freq = write_wavelength_file(wavelength, fname_freq)
    print(f'written {fname_freq}')

    # loop through temperatures to compute hailstone parameters over different
    # melting stages
    prob_break = np.zeros(diam.size)
    h0m = np.nan+np.zeros(diam.size)
    for temp, alt, p_air, dens_air, visc in zip(temp_vec, alt_vec, p_air_vec,
                                                dens_air_vec, visc_vec):
        delta_mw, vel = compute_hail_melting(
            temp, d_ext, d_int, mass_hail, fmw, alpha, delta_h=delta_h,
            dens_air=dens_air, visc=visc, p_air=p_air)
        d_ext, d_int, alpha, fmw, mass_hail, _ = hailstone_change(
            delta_mw, vol_soakable, alpha, fmw, mass_hail)

        h0m, prob_break = compute_breakup_prob(
            h0m, prob_break, fmw, d_ext, alt)

        mass_shed = mass_hail_init - mass_hail

        axis_ratio_int = compute_axis_ratio_dh(d_int)
        axis_ratio_ext = compute_axis_ratio_mh(d_ext, fmw)

        # compute the refractive indices:
        m_water = refractive_index_water(wavelength, temp)

        # internal diameter output can have arbitrary values
        d_int_out = deepcopy(d_int)

        # refractive index of core of melting hailstones
        m_core = np.zeros(diam.size, dtype=complex)
        ind = np.where(fmw == 1.)[0]
        if ind.size > 0:
            # value for completely melted core
            m_core[ind] = m_water
            d_int_out[ind] = 0.001  # set to arbitrary value

        ind = np.where(fmw < 1.)[0]
        if ind.size > 0:
            # values for the core in the process of melting
            m_core[ind] = refractive_index_melting_hail_core(
                wavelength, mass_hail[ind], alpha[ind], fmw[ind],
                d_int_out[ind])

        # refractive index of water coat
        m_water_coat = m_water+np.zeros(diam.size, dtype=complex)
        ind = np.where(fmw == 0.)[0]
        if ind.size > 0:
            # values for a completely frozen hailstone
            m_pure_ice = refractive_index_ice(wavelength, -1)
            m_ice_core = refractive_index_ice_particle(
                m_pure_ice, dens_hail_init)

            m_water_coat[ind] = m_ice_core[ind]
            d_int_out[ind] = 0.001

        # write results
        fname_model_scatt = (
            f'{savedir}sp_{hydro_type}_{args.band}'
            f'_{int(temp*100.):04d}_model_scatt.txt')
        fname_model_scatt = write_melting_hydro_scatt_model(
            d_ext, d_int_out, axis_ratio_ext, axis_ratio_int, m_water_coat,
            m_core, fname_model_scatt)
        print(f'written {fname_model_scatt}')

        fname = (
            f'{savedir}sp_{hydro_type}_{int(temp*100):04d}_model_part.csv')
        fname = write_melting_hail_part_model(
            diam, d_ext, d_int_out, axis_ratio_ext, axis_ratio_int, mass_hail,
            mass_shed, fmw, alpha, prob_break, vel, temp, fname)
        print(f'written {fname}')

        # compute scattering matrix
        if args.compute_scatt:
            fname_scatt = (
                f'{savedir}sp_{hydro_type}_{args.band}'
                f'_{int(temp*100):04d}_tmat.out')
            fname_scatt = compute_tmatrix_2l(
                fname_freq, fname_model_scatt, args.tm2l_dir, fname_scatt)
            print(f'written {fname_scatt}')


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
