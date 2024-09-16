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
from decantation_filtration import DecantationFiltration

def main():

    # Debugging
    make_run   = True
    make_plots = True
    attach_leaching = True

    # Preamble
    end_time = 14.0*unit.day
    time_step = 10.0*unit.minute
    show_time = (True, unit.hour)

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

    if attach_leaching:
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
        plant_net.connect([decant_filt, 'ccd-overflow'], [leaching, 'pre-leach-feed'])
        plant_net.connect([leaching, 'acid-leach-product'], [decant_filt, 'ccd-feed'])
        plant_net.connect([decant_filt, 'std-underflow'], [leaching, 'acid-leach-feed'])

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

        (quant, time_unit) = decant_filt.std_state_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('decant-filt-std-state-liq-volume.png', dpi=300)

        (quant, time_unit) = decant_filt.std_overflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                   #y_label=quant.latex_name+' ['+quant.unit+']')
                   y_label=quant.latex_name+' [kg/min]')
        plt.grid()
        plt.savefig('decant-filt-std-overflow-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.std_underflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                   #y_label=quant.latex_name + ' [' + quant.unit + ']')
                   y_label=quant.latex_name + ' [kg/min]')
        plt.grid()
        plt.savefig('decant-filt-std-underflow-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.ccd_state_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('decant-filt-ccd-state-liq-volume.png', dpi=300)

        (quant, time_unit) = decant_filt.ccd_overflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                   #y_label=quant.latex_name + ' [' + quant.unit + ']')
                   y_label=quant.latex_name + ' [kg/min]')
        plt.grid()
        plt.savefig('decant-filt-ccd-overflow-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.ccd_underflow_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                   #y_label=quant.latex_name + ' [' + quant.unit + ']')
                   y_label=quant.latex_name + ' [kg/min]')
        plt.grid()
        plt.savefig('decant-filt-ccd-underflow-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.filtration_filtrate_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1 / unit.day, x_label='Time [d]',y_scaling = unit.minute,
                   y_label=quant.latex_name + ' [kg/min]')
        plt.grid()
        plt.savefig('decant-filt-filtrate-mass-flowrate.png', dpi=300)

        (quant, time_unit) = decant_filt.filtration_slurry_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1 / unit.day, x_label='Time [d]',y_scaling = unit.minute,
                   y_label=quant.latex_name + ' [kg/min]')
        plt.grid()
        plt.savefig('decant-filt-slurry-mass-flowrate.png', dpi=300)

        if attach_leaching:
        # Leaching plots
            leaching = plant_net.modules[1]

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

            quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                       #y_label=quant.latex_name+' ['+quant.unit+']')
                       y_label=quant.latex_name+' [kg/min]')
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

if __name__ == '__main__':
    main()
