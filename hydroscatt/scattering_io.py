"""
io
==

function for reading/writing

.. autosummary::
    :toctree: generated/

    get_save_dir
    write_melting_hydro_scatt_model
    write_melting_hail_part_model
    write_melting_snow_part_model
    write_wavelength_file
    read_multiple_hydro_part_model
    read_multiple_hydro_scatt_model
    read_melting_hydro_part_model
    read_melting_hydro_scatt_model
    read_scatt_double_layer
    read_wavelength_file

"""

import os
import csv
import glob
from warnings import warn

import numpy as np
import pandas as pd


def get_save_dir(basepath, hydro_type, band, temp, create_dir=True,
                 with_subdirs=False):
    """
    obtains the path to a product directory and eventually creates it

    Parameters
    ----------
    basepath : str
        product base path
    hydro_type : str
        hydrometeor type
    band : str
        frequency band
    temp : float
        temperature
    create_dir : boolean
        If True creates the directory
    with_subdirs : boolean
        If True subdirectories are going to be added to the base path

    Returns
    -------
    savedir : str
        path to product

    """
    if with_subdirs:
        savedir = f'{basepath}{hydro_type}/{band}/t{int(temp*100.):04d}/'
    else:
        savedir = f'{basepath}'

    if create_dir is False:
        return savedir

    if not os.path.isdir(savedir):
        os.makedirs(savedir)

    return savedir


def write_melting_hydro_scatt_model(d_ext, d_int, axis_ratio_ext,
                                    axis_ratio_int, m_ext, m_core, fname):
    """
    write melting hydrometeor model used for 2-layer T-matrix scattering
    computations

    Parameters
    ----------
    d_ext : array of floats
        external diameter (mm)
    d_int : array of floats
        internal core diameter (mm)
    axis_ratio_ext : array of floats
        axis ratio of external shell
    axis_ratio_int : array of floats
        axis ratio of internal core
    m_ext : array of floats
        refractive index of external shell
    m_core : array of floats
        refractive index of inner core
    fname : str
        file name

    Returns
    -------
    fname : str
        name of the file where the data has been written

    """
    d_ext_cm = 0.1*d_ext
    d_int_cm = 0.1*d_int

    err_ext = m_ext.real*m_ext.real-m_ext.imag*m_ext.imag
    eri_ext = 2.*m_ext.real*m_ext.imag

    err_int = m_core.real*m_core.real-m_core.imag*m_core.imag
    eri_int = 2.*m_core.real*m_core.imag

    with open(fname, 'w', newline='', encoding="utf-8") as csvfile:
        csvfile.write(f"{d_ext.size}\n")

        field_names = [
            'd_ext', 'd_int', 'ar_ext', 'ar_int', 'err_ext', 'eri_ext',
            'err_int', 'eri_int']

        writer = csv.DictWriter(csvfile, field_names, delimiter=' ')

        for ind in range(d_ext.size):
            dict_row = {
                'd_ext': d_ext_cm[ind],
                'd_int': d_int_cm[ind],
                'ar_ext': axis_ratio_ext[ind],
                'ar_int': axis_ratio_int[ind],
                'err_ext': err_ext[ind],
                'eri_ext': eri_ext[ind],
                'err_int': err_int[ind],
                'eri_int': eri_int[ind],
            }
            writer.writerow(dict_row)
        csvfile.close()

    return fname


def write_melting_hail_part_model(d_init, d_ext, d_int, axis_ratio_ext,
                                  axis_ratio_int, mass_hs, mass_shed, fmw,
                                  alpha, prob_break, vel, temp, fname):
    """
    write melting hail particle model

    Parameters
    ----------
    d_init : array of floats
        hailstone diameter before the start of the melting process (mm)
    d_ext : array of floats
        external hailstone diameter (mm)
    d_int : array of floats
        internal core hailstone diameter (mm)
    axis_ratio_ext : array of floats
        axis ratio of external part of the hailstone (water coat)
    axis_ratio_int : array of floats
        axis ratio of internal core of the hailstone
    mass_hs : array of floats
        hailstone mass (g)
    mass_shed : array of floats
        mass of shed water (g)
    fmw : array of floats
        fraction mass water of hailstone
    alpha : array of floats or NaNs
        Ratio between the mass of soaked water and the mass of pure ice in
        the hailstone core. The value is NaN if the hailstone is not yet
        soaked
    prob_break : array of floats
        Probability that large drops created by completely melted hailstones
        break at a given altitude
    vel : array of floats
        hailstones terminal velocity
    temp : float
        temperature (deg C)
    fname : str
        file name

    Returns
    -------
    fname : str
        name of the file where the data has been written

    """
    with open(fname, 'w', newline='', encoding="utf-8") as csvfile:
        field_names = [
            'd_init', 'd_ext', 'd_int', 'ar_ext', 'ar_int', 'hail_mass',
            'shed_water_mass', 'fmw', 'alpha', 'prob_break', 'vel']

        csvfile.write(f"{temp}\n")
        writer = csv.DictWriter(csvfile, field_names)
        writer.writeheader()

        for ind in range(mass_hs.size):
            dict_row = {
                'd_init': d_init[ind],
                'd_ext': d_ext[ind],
                'd_int': d_int[ind],
                'ar_ext': axis_ratio_ext[ind],
                'ar_int': axis_ratio_int[ind],
                'hail_mass': mass_hs[ind],
                'shed_water_mass': mass_shed[ind],
                'fmw': fmw[ind],
                'alpha': alpha[ind],
                'prob_break': prob_break[ind],
                'vel': vel[ind]
            }
            writer.writerow(dict_row)
        csvfile.close()

    return fname


def write_melting_snow_part_model(d_init, d_ext, d_int, axis_ratio_ext,
                                  axis_ratio_int, mass, fmw, vel, temp, fname):
    """
    write melting snow particle model

    Parameters
    ----------
    d_init : array of floats
        snowflake diameter before the start of the melting process (mm)
    d_ext : array of floats
        external snowflake diameter (mm)
    d_int : array of floats
        internal core snowflake diameter (mm)
    axis_ratio_ext : array of floats
        axis ratio of external part of the snowflake
    axis_ratio_int : array of floats
        axis ratio of internal core of the snowflake
    mass : array of floats
        snowflake mass (g)
    fmw : array of floats
        fraction mass water of snowflake
    vel : array of floats
        snowflake velocity (m/s)
    temp : float
        temperature (deg C)
    fname : str
        file name

    Returns
    -------
    fname : str
        name of the file where the data has been written

    """
    with open(fname, 'w', newline='', encoding="utf-8") as csvfile:
        field_names = [
            'd_init', 'd_ext', 'd_int', 'ar_ext', 'ar_int', 'mass', 'fmw',
            'vel']

        csvfile.write(f"{temp}\n")
        writer = csv.DictWriter(csvfile, field_names)
        writer.writeheader()

        for ind in range(d_init.size):
            dict_row = {
                'd_init': d_init[ind],
                'd_ext': d_ext[ind],
                'd_int': d_int[ind],
                'ar_ext': axis_ratio_ext[ind],
                'ar_int': axis_ratio_int[ind],
                'mass': mass[ind],
                'fmw': fmw[ind],
                'vel': vel[ind]
            }
            writer.writerow(dict_row)
        csvfile.close()

    return fname


def write_wavelength_file(wavelength, fname):
    """
    write the wavelength file. The wavelength is expressed in cm

    Parameters
    ----------
    wavelength : float
        wavelength (mm)
    fname : str
        file name

    Returns
    -------
    fname : str
        name of the file where the data has been written

    """
    with open(fname, 'w', newline='', encoding="utf-8") as csvfile:
        csvfile.write(f"&ingest\npalamd={wavelength/10.}\n/\n")
        csvfile.close()

    return fname


def read_multiple_hydro_part_model(path, hydro_type, d_max=None):
    """
    read multiple files containing the particle characteristics
    and combines them in a single pandas data frame

    Parameters
    ----------
    path : str
        path to the files containing the particle model data
    hydro_type : str
        hydrometeor type
    d_max : float or None
        if not None the data will be cut at d_max (mm)

    Returns
    -------
    df_model : pandas DataFrame
        df containing the parameters
    temp : float
        the ambient temperature

    """
    flist = glob.glob(f'{path}sp_*{hydro_type}_*_model_part.csv')
    if not flist:
        warn(f'No file founds in {path}sp_*{hydro_type}_*_model_part.csv')
        return None, None

    temp_vec = []
    temp_vec2 = []
    first_file = True
    for fname in flist:
        df_aux, temp = read_melting_hydro_part_model(fname, d_max=d_max)
        temp_vec.append(temp)
        temp_vec2.extend(temp+np.zeros(df_aux.shape[0]))
        if first_file:
            df_model = df_aux
            first_file = False
        else:
            df_model = pd.concat([df_model, df_aux], ignore_index=True)

    temp_vec = np.array(temp_vec)
    temp_vec2 = np.array(temp_vec2)
    df_model['temp'] = temp_vec2
    df_model.sort_values('temp', inplace=True)
    temp_vec = np.sort(temp_vec)

    return df_model, temp_vec


def read_multiple_hydro_scatt_model(path, hydro_type, band='C'):
    """
    read multiple files containing the parameters required to compute
    the scattering properties and combines them in a single pandas data frame

    Parameters
    ----------
    path : str
        path to the files containing the particle model data
    hydro_type : str
        hydrometeor type

    Returns
    -------
    df_model : pandas DataFrame
        df containing the parameters
    temp : float
        the ambient temperature

    """
    flist = glob.glob(f'{path}sp_*{hydro_type}_{band}_*_model_scatt.txt')
    if not flist:
        warn(f'No file founds in {path}sp_*{hydro_type}_{band}_*_'
             f'model_scatt.txt')
        return None, None

    temp_vec = []
    temp_vec2 = []
    first_file = True
    for fname in flist:
        df_aux = read_melting_hydro_scatt_model(fname)
        temp = float(fname.split('_')[-3])/100.
        temp_vec.append(temp)
        temp_vec2.extend(temp+np.zeros(df_aux.shape[0]))
        if first_file:
            df_model = df_aux
            first_file = False
        else:
            df_model = pd.concat([df_model, df_aux], ignore_index=True)

    temp_vec = np.array(temp_vec)
    temp_vec2 = np.array(temp_vec2)
    df_model['temp'] = temp_vec2
    df_model.sort_values('temp', inplace=True)
    temp_vec = np.sort(temp_vec)

    return df_model, temp_vec


def read_melting_hydro_part_model(fname, d_max=None):
    """
    read a file containing parameters describing melting hydrometeors

    Parameters
    ----------
    fname : str
        file name
    d_max : float or None
        if not None the data will be cut at d_max (mm)

    Returns
    -------
    df : pandas DataFrame
        df containing the parameters
    temp : float
        the ambient temperature

    """
    with open(fname, 'r', newline='', encoding="utf-8") as csvfile:
        temp = float(csvfile.readline())
        df = pd.read_csv(csvfile)
        if d_max is not None:
            df = df[df['d_init'] <= d_max]

    return df, temp


def read_melting_hydro_scatt_model(fname, col_names=None, nrows=None):
    """
    read a file containing the parameters required to compute the scattering
    properties of melting hydrometeors

    Parameters
    ----------
    fname : str
        file name
    col_names : list of str or None
        names of each column
    d_max : float or None
        if not None the data will be cut at d_max

    Returns
    -------
    df : pandas DataFrame
        df containing the parameters

    """
    if col_names is None:
        col_names = [
            'd_ext', 'd_int', 'ar_ext', 'ar_int', 'err_ext', 'eri_ext',
            'err_int', 'eri_int']
    if nrows is not None:
        df = pd.read_csv(
            fname, sep=' ', skiprows=1, names=col_names, nrows=nrows)
    else:
        df = pd.read_csv(fname, sep=' ', skiprows=1, names=col_names)
    return df


def read_scatt_double_layer(fname, col_names=None, nrows=None):
    """
    read a file containing the elements of a scattering matrix

    Parameters
    ----------
    fname : str
        file name
    col_names : list of str or None
        names of each column
    nrows : int or None
        number of rows to read

    Returns
    -------
    df : pandas DataFrame
        df containing the parameters

    """
    if col_names is None:
        col_names = [
            'fv180_re', 'fv180_im', 'fh180_re', 'fh180_im', 'fv0_re',
            'fv0_im', 'fh0_re', 'fh0_im']
    if nrows is not None:
        df = pd.read_csv(fname, sep=r'\s+', names=col_names, nrows=nrows)
    else:
        df = pd.read_csv(fname, sep=r'\s+', names=col_names)
    df['fv180'] = -(df['fv180_re']+1j*df['fv180_im'])
    df['fh180'] = df['fh180_re']+1j*df['fh180_im']
    df['fv0'] = df['fv0_re']+1j*df['fv0_im']
    df['fh0'] = df['fh0_re']+1j*df['fh0_im']

    return df


def read_wavelength_file(fname):
    """
    read the wavelength file

    Parameters
    ----------
    fname : str
        file name

    Returns
    -------
    wavelength : float
        wavelength (mm)

    """
    with open(fname, 'r', newline='', encoding="utf-8") as csvfile:
        csvfile.readline()
        wave_text = csvfile.readline()
        csvfile.close()
        wavelength = 10.*float(wave_text.split('=')[1])

    return wavelength
