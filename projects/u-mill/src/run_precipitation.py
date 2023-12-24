#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from precipitation import Precipitation

def main():

    # Debugging
    make_run   = True
    make_plots = True

    # Preamble
    end_time = 1*unit.hour
    time_step = 1.0*unit.minute
    show_time = (True, 5*unit.minute)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Precipitation

    # Steady state conditions for NuSCale case
    #primary_inflow_temp = (320.9+273.15)*unit.kelvin
    #secondary_inflow_temp = (149+273.15)*unit.kelvin

    precipt = Precipitation()  # Create precipitation module

    precipt.name = 'Precipitation'
    precipt.save = True
    precipt.time_step = time_step
    precipt.end_time = end_time
    precipt.show_time = show_time

    plant_net.module(precipt)  # Add precipt module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdow Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Precipitation plots
        precipt = plant_net.modules[0]

        (quant, time_unit) = precipt.precipitation_phase.get_quantity_history('mass-flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipitation-mass-flowrate.png', dpi=300)

        '''
        (quant, time_unit) = precipt.primary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-primary-mass-flowrate.png', dpi=300)

        (quant, time_unit) = precipt.secondary_inflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('precipt-secondary-inflow-temp.png', dpi=300)

        (quant, time_unit) = precipt.secondary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('precipt-secondary-outflow-temp.png', dpi=300)

        (quant, time_unit) = precipt.secondary_inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-secondary-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = precipt.secondary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-secondary-outflow-flowrate.png', dpi=300)

        (quant, time_unit) = precipt.state_phase.get_quantity_history('tau_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-primary-tau.png', dpi=300)

        (quant, time_unit) = precipt.state_phase.get_quantity_history('tau_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-secondary-tau.png', dpi=300)

        (quant, time_unit) = precipt.secondary_outflow_phase.get_quantity_history('quality')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-secondary-quality.png', dpi=300)

        (quant, time_unit) = precipt.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.kilo, x_label='Time [m]',
                   y_label=quant.latex_name+' [k'+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-heatflux.png', dpi=300)

        (quant, time_unit) = precipt.state_phase.get_quantity_history('nusselt_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-nusselt_p.png', dpi=300)

        (quant, time_unit) = precipt.state_phase.get_quantity_history('nusselt_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('precipt-nusselt_s.png', dpi=300)
        '''

if __name__ == '__main__':
    main()
