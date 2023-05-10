"""
graph
=====

Functions to plot scattering parameters

.. autosummary::
    :toctree: generated/

    plot_vertical_profile
    plot_multiple_var
    plot_polvar2
    plot_polvar
    plot_polvar_scat2
    plot_polvar_scat
    plot_psd_scatt_quantities
    plot_sp_scatt_quantities
    plot_psd_scatt_profile
    get_label

"""

from warnings import warn
from matplotlib import pyplot as plt


def plot_vertical_profile(varx_vec, vary, dpi=72,
                          xlabel='external diameter [mm]',
                          ylabel='temperature [C]', titl='plot',
                          fname='./d-temp.png', labels=None,
                          invert_yaxis=True):
    """
    plots the vertical profile of multiple parameters

    Parameters
    ----------
    varx_vec : list of array of floats
        The x variables to plot
    vary : array of floats
        The y variables
    xlabel, ylabel : str
        the x and y axis labels
    titl : str
        plot title
    fname : str
        file name where to save the data
    labels : list of str or None
        if None the plots will have no label. If a list of str these are the
        labels of the plot legend
    invert_yaxis : bool
        if True yaxis will be inverted

    """
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    for ind, varx in enumerate(varx_vec):
        label = None
        if labels is not None:
            label = labels[ind]
        ax.plot(varx, vary, marker='o', linestyle='-', label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)
    if invert_yaxis:
        ax.invert_yaxis()
    ax.grid()
    if labels is not None:
        ax.legend(loc='best')

    fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    print(f'saved {fname}')


def plot_multiple_var(varx_vec, vary_vec, dpi=72,
                      xlabel='external diameter [mm]',
                      ylabel='temperature [C]', titl='plot',
                      fname='./d-temp.png', labels=None, logy=False):
    """
    plots multiple variables in the same plot

    Parameters
    ----------
    varx_vec : list array of floats
        The x variable to plot
    vary_vec : list of array of floats
        The y variables
    xlabel, ylabel : str
        the x and y axis labels
    titl : str
        plot title
    fname : str
        file name where to save the data
    labels : list of str or None
        if None the plots will have no label. If a list of str these are the
        labels of the plot legend

    """
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    for ind, vary in enumerate(vary_vec):
        if len(varx_vec) == len(vary_vec):
            varx = varx_vec[ind]
        else:
            varx = varx_vec[0]
        label = None
        if labels is not None:
            label = labels[ind]
        ax.plot(varx, vary, marker='o', linestyle='-', label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.set_title(titl)
    ax.grid()
    if labels is not None:
        ax.legend(loc='best')

    fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    print(f'saved {fname}')


def plot_polvar2(varx, vary1, vary2, dpi=72, ylabel='reflectivity [dBZ]',
                 xlabel='D0 [mm]', label1='H', label2='V', titl='plot',
                 fname='./dBZ.png'):
    """
    plots 2 variables over the same X axis

    Parameters
    ----------
    varx : array of floats
        The x variable to plot
    vary1, vary2 : array of floats
        The y variables
    ylabel, xlabel : str
        the x and y axis labels
    label1, label2 : str
        the labels of the y variables used in the legend
    titl : str
        plot title
    fname : str
        file name where to save the data

    """
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.plot(varx, vary1, marker='o', linestyle='-', label=label1)
    ax.plot(varx, vary2, marker='o', linestyle='-', label=label2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)
    ax.grid()
    ax.legend(loc='best')

    fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    print(f'saved {fname}')


def plot_polvar(varx, vary, dpi=72,
                xlabel='horizontal reflectivity [dBZ]',
                ylabel='differential reflectivity [dB]', titl='plot',
                fname='./dBZ-ZDR.png'):
    """
    plots 1 variables

    Parameters
    ----------
    varx : array of floats
        The x variable to plot
    vary, vary : array of floats
        The y variable
    ylabel, xlabel : str
        the x and y axis labels
    titl : str
        plot title
    fname : str
        file name where to save the data

    """
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.plot(varx, vary, marker='o', linestyle='-')
    ax.set_title(titl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()

    fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    print(f'saved {fname}')


def plot_polvar_scat2(varx, vary1, vary2, dpi=72, ylabel='reflectivity [dBZ]',
                      xlabel='D0 [mm]', label1='H', label2='V',
                      titl='scatter plot', fname='./dBZ.png'):
    """
    scatter plot of 2 variables over the same X axis

    Parameters
    ----------
    varx : array of floats
        The x variable to plot
    vary1, vary2 : array of floats
        The y variables
    ylabel, xlabel : str
        the x and y axis labels
    label1, label2 : str
        the labels of the y variables used in the legend
    titl : str
        plot title
    fname : str
        file name where to save the data

    """
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.scatter(varx, vary1, marker='o', label=label1)
    ax.scatter(varx, vary2, marker='o', label=label2)
    ax.set_title(titl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.legend(loc='best')

    fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    print(f'saved {fname}')


def plot_polvar_scat(varx, vary, dpi=72, titl='scatter plot',
                     xlabel='horizontal reflectivity [dBZ]',
                     ylabel='differential reflectivity [dB]',
                     fname='./dBZ-ZDR.png'):
    """
    scatter plot of 1 variable

    Parameters
    ----------
    varx : array of floats
        The x variable to plot
    vary : array of floats
        The y variable
    titl : str
        plot title
    ylabel, xlabel : str
        the x and y axis labels
    fname : str
        file name where to save the data

    """
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.scatter(varx, vary, marker='o')
    ax.set_title(titl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()

    fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    print(f'saved {fname}')


def plot_psd_scatt_quantities(df, path, band, temp, hydro_type, ele=0.,
                              x_var_list=['refl_h', 'lwc', 'rr'],
                              y_var_list=['refl', 'ldr', 'A', 'sca_xsect',
                                          'ext_xsect']):
    """
    Plots the selected PSD scattering quantitites

    Parameters
    ----------
    df : Pandas Data Frame
        DataFrame containing the data
    path : str
        path where to store the plots
    band : str
        frequency band
    temp : str or float
        temperature (deg C)
    ele : float
        elevation angle (deg)
    hydro_type : str
        hydrometeor type
    x_var_list, y_var_list : list of str
        list of variables to plot. Must be among:
        D, sca_xsect, ext_xsect, refl, ldr,  zdr, rho_hv, delta_hv, kdp, A, Adp
        If they exist the horizontal and vertical polarized quantity will be
        plot in the same axis

    """
    temp100 = int(temp*100.)
    titl = f'{hydro_type} {band}-band {temp} deg C ele {ele} deg'
    for x_var in x_var_list:
        if x_var not in df:
            warn(f'Unable to plot variable {x_var}')
            continue

        xlabel = get_label(x_var)
        for y_var in y_var_list:
            ylabel = get_label(y_var)
            fname = (
                f'{path}psd_{hydro_type}_{band}_{temp100:04d}'
                f'_ele{int(ele*100.):05d}_{x_var}-{y_var}.png')

            if y_var in ('refl', 'ldr', 'A', 'sca_xsect', 'ext_xsect'):
                if f'{y_var}_h' in df.columns and f'{y_var}_v' in df.columns:
                    if x_var in (f'{y_var}_h', f'{y_var}_v'):
                        if x_var == f'{y_var}_h':
                            y_var_aux = f'{y_var}_v'
                        else:
                            y_var_aux = f'{y_var}_h'
                        ylabel = get_label(y_var_aux)
                        plot_polvar_scat(
                            df[x_var], df[y_var_aux], xlabel=xlabel,
                            ylabel=ylabel, titl=titl, fname=fname)
                    else:
                        plot_polvar_scat2(
                            df[x_var], df[f'{y_var}_h'], df[f'{y_var}_v'],
                            xlabel=xlabel, ylabel=ylabel, label1='H',
                            label2='V', titl=titl, fname=fname)
                else:
                    if f'{y_var}_h' in df.columns:
                        y_var_aux = f'{y_var}_h'
                    elif f'{y_var}_v' in df.columns:
                        y_var_aux = f'{y_var}_v'
                    else:
                        warn(f'Unable to plot variable {y_var}')
                        continue
                    if x_var == y_var_aux:
                        continue

                    ylabel = get_label(y_var_aux)
                    plot_polvar_scat(
                        df[x_var], df[y_var_aux], xlabel=xlabel, ylabel=ylabel,
                        titl=titl, fname=fname)
            else:
                if x_var == y_var:
                    continue
                if y_var not in df:
                    warn(f'Unable to plot variable {y_var}')
                    continue
                plot_polvar_scat(
                    df[x_var], df[y_var], xlabel=xlabel, ylabel=ylabel,
                    titl=titl, fname=fname)


def plot_sp_scatt_quantities(df, path, band, temp, hydro_type, ele=0.,
                             x_var_list=['d'],
                             y_var_list=['refl', 'ldr', 'A', 'sca_xsect',
                                         'ext_xsect']):
    """
    Plots the selected single particle scattering quantitites

    Parameters
    ----------
    df : Pandas Data Frame
        DataFrame containing the data
    path : str
        path where to store the plots
    band : str
        frequency band
    temp : str or float
        temperature (deg C)
    ele : float
        elevation angle (deg)
    hydro_type : str
        hydrometeor type
    x_var_list, y_var_list : list of str
        list of variables to plot. Must be among:
        D, sca_xsect, ext_xsect, refl, ldr,  zdr, rho_hv, delta_hv, kdp, A, Adp
        If they exist the horizontal and vertical polarized quantity will be
        plot in the same axis

    """
    temp100 = int(temp*100.)
    titl = f'{hydro_type} {band}-band {temp} deg C ele {ele} deg'
    for x_var in x_var_list:
        if x_var not in df:
            warn(f'Unable to plot variable {x_var}')
            continue

        xlabel = get_label(x_var)
        for y_var in y_var_list:
            ylabel = get_label(y_var)
            fname = (
                f'{path}sp_{hydro_type}_{band}_{temp100:04d}'
                f'_ele{int(ele*100.):05d}_{x_var}-{y_var}.png')

            if y_var in ('refl', 'ldr', 'A', 'sca_xsect', 'ext_xsect'):
                if f'{y_var}_h' in df.columns and f'{y_var}_v' in df.columns:
                    if x_var in (f'{y_var}_h', f'{y_var}_v'):
                        if x_var == f'{y_var}_h':
                            y_var_aux = f'{y_var}_v'
                        else:
                            y_var_aux = f'{y_var}_h'
                        ylabel = get_label(y_var_aux)
                        plot_polvar(
                            df[x_var], df[y_var_aux], xlabel=xlabel,
                            ylabel=ylabel, titl=titl, fname=fname)
                    else:
                        plot_polvar2(
                            df[x_var], df[f'{y_var}_h'], df[f'{y_var}_v'],
                            xlabel=xlabel, ylabel=ylabel, label1='H',
                            label2='V', titl=titl, fname=fname)
                else:
                    if f'{y_var}_h' in df.columns:
                        y_var_aux = f'{y_var}_h'
                    elif f'{y_var}_v' in df.columns:
                        y_var_aux = f'{y_var}_v'
                    else:
                        warn(f'Unable to plot variable {y_var}')
                        continue

                    ylabel = get_label(y_var_aux)
                    plot_polvar(
                        df[x_var], df[y_var_aux], xlabel=xlabel, ylabel=ylabel,
                        titl=titl, fname=fname)
            else:
                if x_var == y_var:
                    continue
                if y_var not in df:
                    warn(f'Unable to plot variable {y_var}')
                    continue
                plot_polvar(
                    df[x_var], df[y_var], xlabel=xlabel, ylabel=ylabel,
                    titl=titl, fname=fname)


def plot_psd_scatt_profile(df, path, band, hydro_type, ele=0.,
                           x_var_list=['refl', 'ldr', 'A', 'sca_xsect',
                                       'ext_xsect'],
                           y_var_list=['temp']):
    """
    Plots the selected PSD scattering quantitites profile

    Parameters
    ----------
    df : Pandas Data Frame
        DataFrame containing the data
    path : str
        path where to store the plots
    band : str
        frequency band
    ele : float
        elevation angle (deg)
    hydro_type : str
        hydrometeor type
    x_var_list, y_var_list : list of str
        list of variables to plot. Must be among:
        D, sca_xsect, ext_xsect, refl, ldr,  zdr, rho_hv, delta_hv, kdp, A, Adp
        If they exist the horizontal and vertical polarized quantity will be
        plot in the same axis

    """
    titl = f'{hydro_type} {band}-band ele {ele} deg'
    for y_var in y_var_list:
        if y_var not in df:
            warn(f'Unable to plot variable {y_var}')
            continue

        ylabel = get_label(y_var)
        for x_var in x_var_list:
            xlabel = get_label(x_var)
            fname = (
                f'{path}psd_{hydro_type}_{band}_ele{int(ele*100.):05d}'
                f'_{x_var}-{y_var}.png')

            if x_var in ('refl', 'ldr', 'A', 'sca_xsect', 'ext_xsect'):
                if f'{x_var}_h' in df.columns and f'{x_var}_v' in df.columns:
                    plot_vertical_profile(
                        [df[f'{x_var}_h'], df[f'{x_var}_v']], df[y_var],
                        xlabel=xlabel, ylabel=ylabel, labels=['H', 'V'],
                        titl=titl, fname=fname)
                else:
                    if f'{x_var}_h' in df.columns:
                        x_var_aux = f'{x_var}_h'
                    elif f'{x_var}_v' in df.columns:
                        x_var_aux = f'{x_var}_v'
                    else:
                        warn(f'Unable to plot variable {x_var}')
                        continue

                    xlabel = get_label(x_var_aux)
                    plot_vertical_profile(
                        [df[x_var_aux]], df[y_var],
                        xlabel=xlabel, ylabel=ylabel, labels=['H', 'V'],
                        titl=titl, fname=fname)
            else:
                if x_var not in df:
                    warn(f'Unable to plot variable {x_var}')
                    continue
                plot_vertical_profile(
                    [df[x_var]], df[y_var], xlabel=xlabel, ylabel=ylabel,
                    titl=titl, fname=fname)


def get_label(var):
    """
    given a variable identifier returns its name

    Parameters
    ----------
    var : str
        variable identifier

    Returns
    -------
    label : str
       variable name

    """
    if var == 'temp':
        return 'temperature (deg C)'

    if var == 'd':
        return 'equivalent volume diameter (mm)'
    if var == 'd_ext':
        return 'equivalent volume diameter (mm)'
    if var == 'l':
        return 'maximum diameter/length (mm)'

    if var == 'D0':
        return 'D0 (mm)'
    if var == 'Nw':
        return 'Nw'
    if var == 'mu':
        return 'mu'

    if var == 'lwc':
        return '(equivalent) liquid water content (mm3/m3)'
    if var == 'rr':
        return '(equivalent) rainfall rate (mm/h)'

    if var == 'sca_xsect':
        return 'scattering cross-section (dBsm)'
    if var == 'sca_xsect_h':
        return 'scattering cross-section H (dBsm)'
    if var == 'sca_xsect_v':
        return 'scattering cross-section V (dBsm)'

    if var == 'ext_xsect':
        return 'extinction cross-section (dBsm)'
    if var == 'ext_xsect_h':
        return 'extinction cross-section H (dBsm)'
    if var == 'ext_xsect_v':
        return 'extinction cross-section V (dBsm)'

    if var == 'refl':
        return 'reflectivity (dBZ)'
    if var == 'refl_h':
        return 'reflectivity H (dBZ)'
    if var == 'refl_v':
        return 'reflectivity V (dBZ)'

    if var == 'ldr':
        return 'Linear depolarization ratio (dB)'
    if var == 'ldr_h':
        return 'Linear depolarization ratio H (dB)'
    if var == 'ldr_v':
        return 'Linear depolarization ratio V (dB)'

    if var == 'zdr':
        return 'differential reflectivity (dB)'
    if var == 'rho_hv':
        return 'co-polar correlation coefficient'
    if var == 'delta_hv':
        return 'backscattered differential phase (deg)'
    if var == 'kdp':
        return 'specific differential phase (deg/km)'

    if var == 'A':
        return 'specific attenuation (2-way) (dB/km)'
    if var == 'A_h':
        return 'specific attenuation (2-way) H (dB/km)'
    if var == 'A_v':
        return 'specific attenuation (2-way) V (dB/km)'

    if var == 'Adp':
        return 'specific differential attenuation (2-way) (dB/km)'

    warn(f'Unknown label for {var}')
    return ''
