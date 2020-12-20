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

def main():

    # Debugging
    make_plots = True
    make_run   = True

    # Preamble
    end_time = 60*unit.minute
    time_step = 1.5*unit.second
    show_time = (True, 1*unit.minute)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Reactor
    reactor = SMPWR()  # Create reactor module
    reactor.name = "SM-PWR"

    reactor.time_step = time_step
    reactor.end_time = end_time
    reactor.show_time = show_time

    # Steady state condition for NuScale case
    reactor.inflow_cool_temp = unit.convert_temperature(497,'F','K')

    reactor.shutdown = (True, 45*unit.minute)

    plant_net.module(reactor)  # Add reactor module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdow Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Reactor plots
        reactor = plant_net.modules[0]

        (quant, time_unit) = reactor.neutron_phase.get_quantity_history('neutron-dens')
        #quant.plot(x_scaling=1/unit.minute, y_scaling=1/max(quant.value),
        quant.plot(x_scaling=1/unit.minute,
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

        (quant, time_unit) = reactor.state_phase.get_quantity_history('core-temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('reactor-core-temp.png', dpi=300)

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


if __name__ == '__main__':
    main()
