#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from solvex import Solvex

def main():

    # Debugging
    make_run   = True
    make_plots = True

    # Preamble
    end_time = 10.0 * unit.day
    time_step = 10.0 * unit.minute
    show_time = (True, unit.hour)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Solvent extraction

    solvex = Solvex()  # Create solvent extraction module

    solvex.name = 'Solvex'
    solvex.save = True
    solvex.time_step = time_step
    solvex.end_time = end_time
    solvex.show_time = show_time

    plant_net.module(solvex)  # Add solvex module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdow Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Solvent extraction plots
        solvex = plant_net.modules[0]

        (quant, time_unit) = solvex.solvex_state_phase.get_quantity_history('aqueous-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-state-aqueous-volume.png', dpi=300)

        (quant, time_unit) = solvex.solvex_state_phase.get_quantity_history('organic-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-state-organic-volume.png', dpi=300)

        (quant, time_unit) = solvex.solvex_state_phase.get_quantity_history('liquid-volume')

        quant.plot(x_scaling=1/unit.day, x_label='Time [d]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-state-liquid-volume.png', dpi=300)

        (quant, time_unit) = solvex.solvex_raffinate_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.day, y_scaling=unit.minute, x_label='Time [d]',
                   #y_label=quant.latex_name+' ['+quant.unit+']')
                   y_label=quant.latex_name+' [kg/m$3$]')
        plt.grid()
        plt.savefig('solvex-raffinate-mass-flowrate.png', dpi=300)

        '''
        (quant, time_unit) = solvex.primary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('solvex-primary-outflow-temp.png', dpi=300)

        (quant, time_unit) = solvex.primary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-primary-mass-flowrate.png', dpi=300)

        (quant, time_unit) = solvex.secondary_inflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('solvex-secondary-inflow-temp.png', dpi=300)

        (quant, time_unit) = solvex.secondary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('solvex-secondary-outflow-temp.png', dpi=300)

        (quant, time_unit) = solvex.secondary_inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-secondary-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = solvex.secondary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-secondary-outflow-flowrate.png', dpi=300)

        (quant, time_unit) = solvex.state_phase.get_quantity_history('tau_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-primary-tau.png', dpi=300)

        (quant, time_unit) = solvex.state_phase.get_quantity_history('tau_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-secondary-tau.png', dpi=300)

        (quant, time_unit) = solvex.secondary_outflow_phase.get_quantity_history('quality')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-secondary-quality.png', dpi=300)

        (quant, time_unit) = solvex.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.kilo, x_label='Time [m]',
                   y_label=quant.latex_name+' [k'+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-heatflux.png', dpi=300)

        (quant, time_unit) = solvex.state_phase.get_quantity_history('nusselt_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-nusselt_p.png', dpi=300)

        (quant, time_unit) = solvex.state_phase.get_quantity_history('nusselt_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('solvex-nusselt_s.png', dpi=300)
        '''
if __name__ == '__main__':
    main()
