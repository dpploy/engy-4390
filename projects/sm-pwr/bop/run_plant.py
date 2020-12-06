#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import unit
import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

from reactor import SMPWR
from steamer import Steamer
from turbine import Turbine
from condenser import Condenser
from water_heater import WaterHeater

def main():

    # Debugging
    make_plots = True
    make_run   = True

    # Preamble
    end_time = 8*unit.minute
    time_step = 1.5*unit.second
    show_time = (True, 5*unit.minute)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Reactor
    reactor = SMPWR()  # Create reactor module

    reactor.name = 'SMPWR'
    reactor.save = True
    reactor.time_step = time_step
    reactor.end_time = end_time
    reactor.show_time = show_time

    plant_net.module(reactor)  # Add reactor module to network

    # Steamer

    steamer = Steamer()  # Create reactor module

    steamer.name = 'Steamer'
    steamer.save = True
    steamer.time_step = time_step
    steamer.end_time = end_time
    steamer.show_time = show_time

    plant_net.module(steamer)  # Add steamer module to network

    # Turbine

    turbine = Turbine()  # Create reactor module

    turbine.name = 'Turbine'
    turbine.save = True
    turbine.time_step = time_step
    turbine.end_time = end_time
    turbine.show_time = show_time

    plant_net.module(turbine)  # Add steamer module to network

    '''Condenser'''

    condenser = Condenser()  # Create condenser module

    condenser.name = 'Condenser'
    condenser.save = True
    condenser.time_step = time_step
    condenser.end_time = end_time
    condenser.show_time = show_time

    plant_net.module(condenser)  # Add condenser module to network`

    '''Feedwater Heating system'''

    water_heater = WaterHeater()  # Create water_heater module

    water_heater.name = 'Water Heater'
    water_heater.save = True
    water_heater.time_step = time_step
    water_heater.end_time = end_time
    water_heater.show_time = show_time

    plant_net.module(water_heater)  # Add water_heater module to network

    # Balance of Plant Network Connectivity

    plant_net.connect([reactor, 'coolant-outflow'], [steamer, 'primary-inflow'])
    plant_net.connect([steamer, 'primary-outflow'], [reactor, 'coolant-inflow'])
    plant_net.connect([steamer, 'secondary-outflow'], [turbine, 'inflow'])
    plant_net.connect([turbine, 'outflow'], [condenser, 'inflow'])
    #plant_net.connect([turbine, 'process-heat'], [water_heater, 'heat'])
    plant_net.connect([condenser, 'outflow'], [water_heater, 'inflow'])
    plant_net.connect([water_heater, 'outflow'], [steamer, 'secondary-inflow'])

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Reactor plots
        reactor = plant_net.modules[0]

        '''
        (quant, time_unit) = reactor.neutron_phase.get_quantity_history('neutron-dens')
        quant.plot(x_scaling=1/unit.minute, y_scaling=1/max(quant.value),
                   x_label='Time [m]', y_label=quant.latex_name+' ['+quant.unit+']')

        plt.grid()
        plt.savefig('reactor-neutron-dens.png', dpi=300)

        (quant, time_unit) = reactor.neutron_phase.get_quantity_history('delayed-neutrons-cc')
        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')

        plt.grid()
        plt.savefig('reactor-delayed-neutrons-cc.png', dpi=300)

        (quant, time_unit) = reactor.coolant_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')

        plt.grid()
        plt.savefig('reactor-coolant-outflow-temp.png', dpi=300)

        (quant, time_unit) = reactor.reactor_phase.get_quantity_history('fuel-temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('reactor-fuel-temp.png', dpi=300)
        '''

        (quant, time_unit) = reactor.state_phase.get_quantity_history('power')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-power.png', dpi=300)

        (quant, time_unit) = reactor.state_phase.get_quantity_history('reynolds')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+r' [$\times 10^6$'+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-reynolds.png', dpi=300)

        (quant, time_unit) = reactor.state_phase.get_quantity_history('prandtl')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+r' ['+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-prandtl.png', dpi=300)

        (quant, time_unit) = reactor.coolant_outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+r' ['+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-mass-flowrate.png', dpi=300)

        (quant, time_unit) = reactor.state_phase.get_quantity_history('heatflux')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-heatflux.png', dpi=300)

        (quant, time_unit) = reactor.state_phase.get_quantity_history('inlet-temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('reactor-coolant-inflow-temp.png', dpi=300)

        (quant, time_unit) = reactor.state_phase.get_quantity_history('nusselt')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-nusselt.png', dpi=300)
        '''
        (quant, time_unit) = reactor.state_phase.get_quantity_history('tau')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-coolant-tau.png', dpi=300)

        (quant, time_unit) = reactor.coolant_outflow_phase.get_quantity_history('quality')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+r' ['+quant.unit+']')
        plt.grid()
        plt.savefig('reactor-coolant-outflow-quality.png', dpi=300)
        '''

        # Steamer plots
        steamer = plant_net.modules[1]

        (quant, time_unit) = steamer.primary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-primary-outflow-temp.png', dpi=300)

        (quant, time_unit) = steamer.secondary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-secondary-outflow-temp.png', dpi=300)

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
        plt.savefig('steamer-nusselt_s.png', dpi=300)


        # Turbine plots
        turbine = plant_net.modules[2]

        (quant, time_unit) = turbine.state_phase.get_quantity_history('power')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('turbine-power.png', dpi=300)

        # Condenser plots
        condenser = plant_net.modules[3]

        (quant, time_unit) = condenser.outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('condenser-outflow-temp.png', dpi=300)

        # Water heater plots
        water_heater = plant_net.modules[4]

        (quant, time_unit) = water_heater.outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('water_heater-outflow-temp.png', dpi=300)

    plant.close()  # Properly shutdow plant

if __name__ == '__main__':
    main()
