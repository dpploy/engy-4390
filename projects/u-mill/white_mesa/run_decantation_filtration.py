#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from leaching import Leaching
from decantation_filtration import DecantationFiltration

def main():

    # Debugging
    make_run   = True
    make_plots = True

    # Preamble
    end_time = 7.0*unit.hour
    time_step = 1.0*unit.minute
    show_time = (True, 5*unit.minute)

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

    # Leaching

    leaching = Leaching() # Create a Leaching module

    leaching.name = 'Leaching'
    leaching.save = True
    leaching.time_step = time_step
    leaching.end_time = end_time
    leaching.show_time = show_time

    plant_net.module(leaching)

    # Balance of Plant Network Connectivity

    plant_net.connect([leaching, 'pre-leach-product'], [decant_filt, 'std-feed'])

    #plant_net.draw(engine='circo', node_shape='folder')
    plant_net.draw(engine='dot', node_shape='folder', size='600,1200')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdown Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Decantation plots
        decant_filt = plant_net.modules[0]

        (quant, time_unit) = decant_filt.single_tank_decantation_overflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('decant-filt-std-overflow-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.single_tank_decantation_underflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1 / unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name + ' [' + quant.unit + ']')
        plt.grid()
        plt.savefig('decant-filt-std-underflow-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.ccd_underflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1 / unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name + ' [' + quant.unit + ']')
        plt.grid()
        plt.savefig('decant-filt-ccd-underflow-mass-flowrate.png', dpi=300)

        # Leaching plots
        leaching = plant_net.modules[1]

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-product-mass-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('mass-density')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-product-mass-density.png', dpi=300)

if __name__ == '__main__':
    main()
