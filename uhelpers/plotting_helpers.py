"""Helper functions for recurring plotting tasks

Authors
-------

    Johannes Sahlmann

Use
---

"""
import os
import numpy as np
import pylab as pl
from scipy.stats import norm


def histogram_with_gaussian_fit(omc, facecolors=None, labels=None, titles=None, linecolors=None, xlabel='value', normed=0,
                                save_plot=0, out_dir='', name_seed='', separate_panels=False, show_fit=True, **kwargs):
    """Plot one or several histograms and perform Gaussian fit(s).

    Parameters
    ----------
    omc
    facecolors
    labels
    titles
    linecolors
    xlabel
    normed
    save_plot
    out_dir
    name_seed
    separate_panels
    show_fit
    kwargs

    Returns
    -------

    """
    if omc.ndim == 1:
        Nhist = 1
        omc = np.expand_dims(omc, axis=-1)
    else:
        Nhist = omc.shape[1]

    if facecolors is None:
        facecolors = ['grey'] * Nhist
        linecolors = ['k'] * Nhist
    if labels is None:
        labels = ['data %d' % j for j in np.arange(Nhist)]

    if separate_panels:
        fig = pl.figure(figsize=(12, 5), facecolor='w', edgecolor='k');
        pl.clf()
        alpha = 0.8
        linecolors = ['k'] * Nhist
    else:
        fig = pl.figure(figsize=(8, 6), facecolor='w', edgecolor='k');
        pl.clf()
        alpha = 0.5

    for i in np.arange(Nhist):
        if separate_panels:
            pl.subplot(1, 2, i + 1)

        data = omc[:, i]
        # from http://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python
        (mu, sigma) = norm.fit(data)
        if show_fit == False:
            histlabel = labels[i]
        else:
            histlabel = None
        n, bins, patches = pl.hist(data, normed=normed, facecolor=facecolors[i], color=linecolors[i], alpha=alpha,
                                   histtype='stepfilled', label=histlabel, **kwargs)
        if normed:
            normFact = 1.
            ylabel = 'Probability'
        else:
            normFact = np.sum(n) * np.mean(np.diff(bins));
            ylabel = 'N'
        if show_fit:
            y = norm.pdf(bins, mu, sigma)  # add a 'best fit' line
            l = pl.plot(bins, y * normFact, 'k-', linewidth=2, color=linecolors[i],
                        label='{0:s}: $\mu$={2:1.3f}$\pm${1:1.3f}'.format(labels[i], sigma, mu))

        if titles is not None:
            if type(titles) == list:
                pl.title(titles[i])
            else:
                pl.title(titles)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)  # pl.ylim((0,max(n)+1))
        pl.legend(loc='best')
        pl.show()

    fig.tight_layout(h_pad=0.0)
    if save_plot:
        figName = os.path.join(out_dir, '%s_distortionFit_residualWithFit.pdf' % name_seed)
        pl.savefig(figName, transparent=True, bbox_inches='tight', pad_inches=0)


def multiple_histograms(all_data, facecolors=None, labels=None, titles=None,
                             linecolors=None, xlabel='value', normed=0, save_plot=0, out_dir='',
                             name_seed='', separate_panels=False, show_fit=True, **kwargs):
    """

    Parameters
    ----------
    all_data
    facecolors
    labels
    titles
    linecolors
    xlabel
    normed
    save_plot
    out_dir
    name_seed
    separate_panels
    show_fit
    kwargs

    Returns
    -------

    """
    Nhist = len(all_data)

    if facecolors is None:
        facecolors = ['grey'] * Nhist
    if linecolors is None:
        linecolors = ['k'] * Nhist
    if labels is None:
        labels = ['data %d' % j for j in np.arange(Nhist)]

    if separate_panels:
        fig = pl.figure(figsize=(12, 5), facecolor='w', edgecolor='k')
        pl.clf()
        alpha = 0.8
        linecolors = ['k'] * Nhist
    else:
        fig = pl.figure(figsize=(8, 4), facecolor='w', edgecolor='k')
        pl.clf()
        alpha = 0.5

    for i in np.arange(Nhist):
        if separate_panels:
            pl.subplot(1, 2, i + 1)

        data = all_data[i]
        histlabel = labels[i]
        n, bins, patches = pl.hist(data, normed=normed, facecolor=facecolors[i],
                                   color=linecolors[i], alpha=alpha, histtype='stepfilled',
                                   label=histlabel, **kwargs)
        if normed:
            normFact = 1.
            ylabel = 'Probability'
        else:
            normFact = np.sum(n) * np.mean(np.diff(bins))
            ylabel = 'N'

        if titles is not None:
            pl.title(titles[i])
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)  # pl.ylim((0,max(n)+1))
        pl.legend(loc=2)
        pl.show()

    fig.tight_layout(h_pad=0.0)
