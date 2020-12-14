#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Balance of plant helper functions.
"""

from matplotlib import pyplot as plt

def plot(dict, n_cols=3):
    """Auxiliar help function to plot phase quantities in a Cortix module.

    The particular phase in the Cortix module that has an `actor` matching an
    element of `quant_names` will be searched for. If found, a plot will be created.

    Parameters
    ----------
    module: Module
        A Cortix Module object.
    quant_names: list(str)
        List of names of phase variables to plot.

    """

    plt.rcParams['figure.figsize'] = [20, 4]
    plt.figure(1)

    n_plots = 0
    for (phase,qnames) in dict.items():
        n_plots += len(qnames)

    if n_plots > n_cols:
        n_rows = int((n_plots-(n_plots%n_cols))/n_cols)
    else:
        n_rows = 1

    # Find phase variable names in module
    #phases = list()
    #for i in dir(module):
    #    if i[0] != '_' and 'phase' in i:
    #        phases.append(i)

    # Get the data stored in the Cortix phases
    #iplot = 1
    #for qname in quant_names:
    #    if qname in eval(phase+'.actors'):
    #            (quant,_) = eval(module+'.'+phase+'.get_quantity_history("'+qname+'")')
    #            plt.subplot(nrows,ncols,iplot)
    #            plt.plot(quant.value.index, quant.value.values)
    #            plt.title(quant.info)
    #            plt.grid(True)
    #            iplot += 1

    # Get the data stored in the Cortix phases
    iplot = 1
    for (phase, qnames) in dict.items():
        for qname in qnames:
            if qname in eval(phase+'.actors'):
                (quant,_) = eval(phase+'.get_quantity_history("'+qname+'")')
                plt.subplot(n_rows, n_cols, iplot)
                plt.plot(quant.value.index, quant.value.values)
                plt.title(quant.info)
                plt.grid(True)
                iplot += 1
