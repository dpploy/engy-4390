#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from decantation_filtration import DecantationFiltration

def main():

    # Debugging
    make_run   = True
    make_plots = False

    # Preamble
    end_time = 1.0*unit.minute
    time_step = 1.0*unit.second
    show_time = (True, 20*unit.second)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Decantation-filtration

    decant_filt = DecantationFiltration()  # Create decantation/filtration module
    decant_filt.name = 'Decantation-Filtration'
    decant_filt.save = True
    decant_filt.time_step = time_step
    decant_filt.end_time = end_time
    decant_filt.show_time = show_time

    plant_net.module(decant_filt)  # Add filtration module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdown Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # filtration plots
        filtration = plant_net.modules[0]

        (quant, time_unit) = filtration.single_tank_decantation_raffinate_feed_solids_massfrac.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('filtration-primary-outflow-temp.png', dpi=300)

        '''
        (quant, time_unit) = filtration.primary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1 / unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name + ' [C]')
        plt.grid()
        plt.savefig('filtration-primary-outflow-temp.png', dpi=300)
        '''

if __name__ == '__main__':
    main()
