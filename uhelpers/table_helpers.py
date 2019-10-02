"""Helper functions that work on astropy tables

Authors
-------

    Johannes Sahlmann

Use
---

"""

import os
import numpy as np
import pylab as pl

def print_column_statistics(t):
    """
    print basic statistics of columns in astropy table
    """
    for colnam in t.colnames:
        try:
            print('Mean %s : %3.3f   median %3.3f   std %3.3f' % (colnam, np.mean(t[colnam]), np.median(t[colnam]), np.std(t[colnam])))
        except:
            pass


def plot_columns_simple(t, plot_dir, save_plot=True, name_seed='table', selected_columns=None,
                        overwrite=False, show_plot=False, highlight_index=None, units=None):
    """Plot astropy column values versus row number.

    :param t:
    :param plot_dir:
    :return:
    """
    for colnam in t.colnames:
        fig_name = os.path.join(plot_dir, '%s_%s.pdf' % (name_seed, colnam))
        if (not os.path.isfile(fig_name)) | (overwrite) | (show_plot):

            if selected_columns is None:
                pass
            elif (selected_columns is not None) & (colnam not in selected_columns):
                continue
            try:
                if t[colnam].dtype == 'object':
                    coord = np.array(t[colnam]).astype(np.float)
                else:
                    coord = np.array(t[colnam])

                if units is None:
                    unit = t[colnam].unit
                else:
                    unit = units[colnam]
                pl.figure(figsize=(6, 6), facecolor='w', edgecolor='k'); pl.clf();
                pl.plot(coord, 'bo')
                if highlight_index is not None:
                    pl.plot(highlight_index, coord[highlight_index], 'go', ms=20)
                pl.title('%s (%s)' % (colnam, unit))
                pl.ylabel('%s (%s)' % (colnam, unit))
                if show_plot:
                    pl.show()
                if save_plot:
                    if not os.path.isdir(plot_dir):
                        os.makedirs(plot_dir)
                    fig_name = os.path.join(plot_dir, '%s_%s.pdf' % (name_seed, colnam))
                    pl.savefig(fig_name, transparent=True, bbox_inches='tight', pad_inches=0.05, overwrite=overwrite)
            except:
                print('Exception occurred')
                pass
