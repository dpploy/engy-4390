#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from evaporation_calcination import EvaporationCalcination

def main():

    # Debugging
    make_run   = True
    make_plots = True

    # Preamble
    end_time = 1.0*unit.minute
    time_step = 1.0*unit.second
    show_time = (True, 20*unit.second)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Evaporation-calcination

    evap_calc = EvaporationCalcination()  # Create reactor module

    evap_calc.name = 'Evaporation-Calcination'
    evap_calc.save = True
    evap_calc.time_step = time_step
    evap_calc.end_time = end_time
    evap_calc.show_time = show_time

    plant_net.module(evap_calc)  # Add steamer module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdow Cortix

    # Evaporation and Calcination Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        evaporation_calcination = plant_net.modules[0]

        (quant, time_unit) = evaporation_calcination.calcination_product_phase.get_quantity_history('mass_flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('Calcination-u3o8-product-mass-flowrate.png', dpi=300)

        '''(quant, time_unit) = leaching.acidleach_phase.get_quantity_history('mass_flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('leaching-acidleach-mass-flowrate.png', dpi=300)'''
        
        
        '''steamer = plant_net.modules[0]

        (quant, time_unit) = steamer.primary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-primary-outflow-temp.png', dpi=300)

        (quant, time_unit) = steamer.primary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-primary-mass-flowrate.png', dpi=300)

        (quant, time_unit) = steamer.secondary_inflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-secondary-inflow-temp.png', dpi=300)

        (quant, time_unit) = steamer.secondary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-secondary-outflow-temp.png', dpi=300)

        (quant, time_unit) = steamer.secondary_inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-secondary-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = steamer.secondary_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-secondary-outflow-flowrate.png', dpi=300)

        (quant, time_unit) = steamer.state_phase.get_quantity_history('tau_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-primary-tau.png', dpi=300)

        (quant, time_unit) = steamer.state_phase.get_quantity_history('tau_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-secondary-tau.png', dpi=300)

        (quant, time_unit) = steamer.secondary_outflow_phase.get_quantity_history('quality')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-secondary-quality.png', dpi=300)

        (quant, time_unit) = steamer.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.kilo, x_label='Time [m]',
                   y_label=quant.latex_name+' [k'+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-heatflux.png', dpi=300)

        (quant, time_unit) = steamer.state_phase.get_quantity_history('nusselt_p')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-nusselt_p.png', dpi=300)

        (quant, time_unit) = steamer.state_phase.get_quantity_history('nusselt_s')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-nusselt_s.png', dpi=300)'''

if __name__ == '__main__':
    main()
