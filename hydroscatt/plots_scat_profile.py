#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
plots_scat_profile
================================================

Auxiliary vertical profile plots of scattering properties

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import pandas as pd

from scattering_io import get_save_dir
from graph import plot_psd_scatt_profile


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
        '--fname', type=str,
        default='psd_profile_melting_hail_C_ele00000_scattering.csv',
        help='name of file containing the scattering parameters')

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

    psd_x_var_list = [
        'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']
    psd_y_var_list = ['temp']

    savedir = get_save_dir(
        args.path, None, None, None, create_dir=True,
        with_subdirs=False)

    df_psd = pd.read_csv(f'{savedir}{args.fname}')
    plot_psd_scatt_profile(
        df_psd, savedir, args.band, args.hydro_type,
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
