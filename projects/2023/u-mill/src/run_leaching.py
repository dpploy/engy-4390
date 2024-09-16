#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network
from cortix import Units as unit

from leaching import Leaching

def main():

    # Debugging
    make_run   = True
    make_plots = True

    # Preamble
    end_time = 4*unit.day
    time_step = 10.0*unit.minute
    show_time = (True, unit.hour)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Leaching

    leaching = Leaching()  # Create reactor module

    leaching.name = 'Leaching'
    leaching.save = True
    leaching.time_step = time_step
    leaching.end_time = end_time
    leaching.show_time = show_time

    plant_net.module(leaching)  # Add leaching module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder', ports=True)

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdown Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Leaching plots
        leaching = plant_net.modules[0]

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                   #y_label=quant.latex_name+' ['+quant.unit+']')
                   y_label=quant.latex_name+' [kg/min]')
        plt.grid()
        plt.savefig('leaching-preleach-mass-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('mass-density')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-mass-density.png', dpi=300)

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-liq-volume.png', dpi=300)

        (quant, time_unit) = leaching.acidleach_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-mass-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.acidleach_phase.get_quantity_history('mass-density')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-mass-density.png', dpi=300)

        (quant, time_unit) = leaching.acidleach_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-liq-volume.png', dpi=300)

        '''
        (quant, time_unit) = leaching.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.kilo, x_label='Time [m]',
                   y_label=quant.latex_name+' [k'+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-heatflux.png', dpi=300)
        '''

if __name__ == '__main__':
    main()
