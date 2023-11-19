#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

import unit

from cortix import Cortix
from cortix import Network

from clarification_filtration import ClarificationFiltration #Will have to add a filtration class in filtration.py file

def main():

    # Debugging
    make_plots = False
    make_run   = False

    # Preamble
    end_time = 10*unit.minute
    time_step = 1.5*unit.second
    show_time = (True, 5*unit.minute)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Filtration

    cfiltration = ClarificationFiltration()  # Create Clarification/filtration module

    # Steady state conditions for NuSCale case
    #primary_inflow_temp = (320.9+273.15)*unit.kelvin
    #secondary_inflow_temp = (149+273.15)*unit.kelvin
    #filtration = filtration(primary_inflow_temp, secondary_inflow_temp)  # Create reactor module
    '''
    cfiltration.name = 'filtration'
    cfiltration.save = True
    cfiltration.time_step = time_step
    cfiltration.end_time = end_time
    cfiltration.show_time = show_time

    plant_net.module(filtration)  # Add filtration module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')
    '''
    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdown Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # filtration plots
        filtration = plant_net.modules[0]

        (quant, time_unit) = filtration.primary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('filtration-primary-outflow-temp.png', dpi=300)

        (quant, time_unit) = filtration.primary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-primary-mass-flowrate.png', dpi=300)

        (quant, time_unit) = filtration.secondary_inflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('filtration-secondary-inflow-temp.png', dpi=300)

        (quant, time_unit) = filtration.secondary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('filtration-secondary-outflow-temp.png', dpi=300)

        (quant, time_unit) = filtration.secondary_inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-secondary-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = filtration.secondary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-secondary-outflow-flowrate.png', dpi=300)

        (quant, time_unit) = filtration.state_phase.get_quantity_history('tau_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-primary-tau.png', dpi=300)

        (quant, time_unit) = filtration.state_phase.get_quantity_history('tau_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-secondary-tau.png', dpi=300)

        (quant, time_unit) = filtration.secondary_outflow_phase.get_quantity_history('quality')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-secondary-quality.png', dpi=300)

        (quant, time_unit) = filtration.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.kilo, x_label='Time [m]',
                   y_label=quant.latex_name+' [k'+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-heatflux.png', dpi=300)

        (quant, time_unit) = filtration.state_phase.get_quantity_history('nusselt_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-nusselt_p.png', dpi=300)

        (quant, time_unit) = filtration.state_phase.get_quantity_history('nusselt_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('filtration-nusselt_s.png', dpi=300)

if __name__ == '__main__':
    main()
