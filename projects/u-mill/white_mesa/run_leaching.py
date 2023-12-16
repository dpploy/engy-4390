#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from leaching import Leaching #Will have to add a leaching class in leaching.py file

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

    # Leaching

    leaching = Leaching()  # Create reactor module

    leaching.name = 'Leaching'
    leaching.save = True
    leaching.time_step = time_step
    leaching.end_time = end_time
    leaching.show_time = show_time

    plant_net.module(leaching)  # Add leaching module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdown Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Leaching plots
        leaching = plant_net.modules[0]

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.hour, x_label='Time [h]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-mass-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('mass-density')

        quant.plot(x_scaling=1/unit.hour, x_label='Time [h]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-mass-density.png', dpi=300)

        (quant, time_unit) = leaching.preleach_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.hour, x_label='Time [h]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-preleach-liq-volume.png', dpi=300)

        (quant, time_unit) = leaching.acidleach_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.hour, x_label='Time [h]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-mass-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.acidleach_phase.get_quantity_history('mass-density')

        quant.plot(x_scaling=1/unit.hour, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-mass-density.png', dpi=300)

        (quant, time_unit) = leaching.acidleach_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.hour, x_label='Time [h]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-liq-volume.png', dpi=300)
        '''

        (quant, time_unit) = leaching.secondary_inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-secondary-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.secondary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-secondary-outflow-flowrate.png', dpi=300)

        (quant, time_unit) = leaching.state_phase.get_quantity_history('tau_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-primary-tau.png', dpi=300)

        (quant, time_unit) = leaching.state_phase.get_quantity_history('tau_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-secondary-tau.png', dpi=300)

        (quant, time_unit) = leaching.secondary_outflow_phase.get_quantity_history('quality')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-secondary-quality.png', dpi=300)

        (quant, time_unit) = leaching.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.kilo, x_label='Time [m]',
                   y_label=quant.latex_name+' [k'+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-heatflux.png', dpi=300)

        (quant, time_unit) = leaching.state_phase.get_quantity_history('nusselt_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-nusselt_p.png', dpi=300)

        (quant, time_unit) = leaching.state_phase.get_quantity_history('nusselt_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-nusselt_s.png', dpi=300)
        '''

if __name__ == '__main__':
    main()
