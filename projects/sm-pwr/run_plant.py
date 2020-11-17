#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import scipy.constants as unit

from cortix import Cortix
from cortix import Network

from reactor import SMPWR
from steamer import Steamer

def main():

    end_time = 1 * unit.hour
    unit.second = 1.0
    time_step = 30.0 * unit.second
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

    # Balance of Plant Network Connectivity

    plant_net.connect([reactor, 'coolant-outflow'], [steamer, 'primary-inflow'])
    plant_net.connect([steamer, 'primary-outflow'], [reactor, 'coolant-inflow'])

    plant_net.draw()

    # Run

    plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdow plant

    #reactor = plant_net.modules[0]

    #(quant, time_unit) = reactor.neutron_phase.get_quantity_history('neutron-dens')
    #quant.plot(x_scaling=1/unit.minute, y_scaling=1/max(quant.value), x_label='Time [m]', y_label=quant.formal_name+' ['+quant.unit+']')

    #(quant, time_unit) = reactor.neutron_phase.get_quantity_history('delayed-neutrons-cc')
    #quant.plot(x_scaling=1/unit.minute, x_label='Time [m]', y_label=quant.formal_name+' ['+quant.unit+']')

    #(quant, time_unit) = reactor.reactor_phase.get_quantity_history('fuel-temp')
    #quant.plot(x_scaling=1/unit.minute, x_label='Time [m]', y_label=quant.formal_name+' ['+quant.unit+']')

    #(quant, time_unit) = reactor.coolant_outflow_phase.get_quantity_history('temp')
    #quant.plot(x_scaling=1/unit.minute, x_label='Time [m]', y_label=quant.formal_name+' ['+quant.unit+']')

if __name__ == '__main__':
    main()
