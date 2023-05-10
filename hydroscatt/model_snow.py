#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
model_snow
================================================

modelling of the characteristics of snow

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
from copy import deepcopy

import numpy as np

from pytmatrix import tmatrix_aux

from part_descrip import compute_axis_ratio_snow, compute_snow_mass
from part_descrip import compute_snow_density, compute_snow_melting
from part_descrip import compute_equi_vol_diam, snowflake_change

from refractivity import refractive_index_water
from refractivity import refractive_index_ice
from refractivity import refractive_index_melting_snow_core
from refractivity import refractive_index_melting_snow_shell

from scattering_io import get_save_dir, write_melting_hydro_scatt_model
from scattering_io import write_melting_snow_part_model, write_wavelength_file

from scattering import compute_tmatrix_2l

print(__doc__)


def main():
    """
    main
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to snow modelling framework')

    # keyword arguments
    parser.add_argument(
        '--path', type=str,
        default='/utemp/mdso/figuerasiventuraj/hydroscatt_products/',
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

    args = parser.parse_args()

    print(f'====== snow modelling started: '
          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
    atexit.register(_print_end_msg,
                    "====== snow modelling finished: ")

    with_subdirs = False
    write_intermediate_results = True

    hydro_type = 'melting_snow'

    delta_h = 10.  # m
    if args.temp > 0:
        temp_min = 0.
        lapse_rate = 6.5  # deg C/km
        step = delta_h*lapse_rate/1000.
        temp_vec = np.arange(temp_min, args.temp+step, step)
        temp_ice = -1
    else:
        temp_vec = np.array([args.temp])
        temp_ice = args.temp

    d_max = 20.
    d_min = 0.1
    step = 0.1
    d_ds = np.arange(d_min, d_max+step, step)

    if args.band == 'S':
        wavelength = tmatrix_aux.wl_S
    elif args.band == 'C':
        wavelength = tmatrix_aux.wl_C
    elif args.band == 'X':
        wavelength = tmatrix_aux.wl_X

    print(hydro_type)
    print('band:', args.band)
    print('temp:', args.temp)

    mass_snow = compute_snow_mass(d_ds)
    dens_ds, dens_ds_core, dens_ds_shell = compute_snow_density(
        d_ds, mass_snow, alpha=0.5)

    d_rd = compute_equi_vol_diam(mass_snow)  # equivalent raindrop diameter
    ar_rain = np.zeros(d_rd.size)
    for ind, d_drop in enumerate(d_rd):
        ar_rain[ind] = tmatrix_aux.dsr_thurai_2007(d_drop)
    ar_ds = compute_axis_ratio_snow(d_ds)

    d_ms = deepcopy(d_ds)
    ar_ms = deepcopy(ar_ds)
    fmw = np.zeros(d_ds.size)

    savedir = get_save_dir(
        args.path, hydro_type, args.band, args.temp, create_dir=True,
        with_subdirs=with_subdirs)

    # loop through temperatures to compute hailstone parameters over different
    # melting stages
    for temp in temp_vec:
        delta_mw, vel = compute_snow_melting(
            temp, d_ds, d_ms, mass_snow, ar_ms, fmw, delta_h=delta_h)

        (d_ms, d_core, ar_ms, _, dens_core, dens_shell,
         fmw) = snowflake_change(
            delta_mw, mass_snow, d_ds, dens_ds, dens_ds_core, dens_ds_shell,
            fmw, ar_ds, ar_rain)

        if write_intermediate_results:
            fname = (
                f'{args.path}sp_{hydro_type}_{int(temp*100):04d}'
                f'_model_part.csv')
            fname = write_melting_snow_part_model(
                d_ds, d_ms, d_core, ar_ms, ar_ds, mass_snow, fmw, vel, temp,
                fname)
            print(f'written {fname}')

    # compute refractive index
    m_pure_ice = refractive_index_ice(wavelength, temp_ice)
    m_water = refractive_index_water(wavelength, 0.)
    m_air = 1.+1j*0.

    m_core = refractive_index_melting_snow_core(
        m_pure_ice, m_water, m_air, dens_ds_core, dens_core, fmw, d_core)
    m_shell = refractive_index_melting_snow_shell(
        m_pure_ice, m_water, m_air, dens_ds_shell, dens_shell, fmw, d_core,
        d_ms)

    ind = np.where(fmw == 1.)[0]
    if ind.size > 0:
        d_core[ind] = 0.001
        m_core[ind] = m_water
        m_shell[ind] = m_water

    # write results
    fname_model_scatt = (
        f'{savedir}sp_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
        f'_model_scatt.txt')
    fname_model_scatt = write_melting_hydro_scatt_model(
        d_ms, d_core, ar_ms, ar_ds, m_shell, m_core, fname_model_scatt)
    print(f'written {fname_model_scatt}')

    fname_freq = (
        f'{savedir}sp_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
        f'_freq.inp')
    fname_freq = write_wavelength_file(wavelength, fname_freq)
    print(f'written {fname_freq}')

    fname = (
        f'{savedir}sp_{hydro_type}_{int(args.temp*100):04d}_model_part.csv')
    fname = write_melting_snow_part_model(
        d_ds, d_ms, d_core, ar_ms, ar_ds, mass_snow, fmw, vel, args.temp,
        fname)
    print(f'written {fname}')

    # compute scattering matrix
    if args.compute_scatt:
        fname_scatt = (
            f'{savedir}sp_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
            f'_tmat.out')
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
