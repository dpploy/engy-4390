#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network

import unit

from leaching import Leaching
from decantation_filtration import DecantationFiltration
from solvex import Solvex
from precipitation import Precipitation
from evaporation_calcination import EvaporationCalcination

def main():

    # Debugging
    make_run   = False
    make_plots = False

    # Preamble
    end_time = 10*unit.minute
    time_step = 1*unit.second
    show_time = (True, 1*unit.minute)

    whte_mesa = Cortix(use_mpi=False, splash=True) # System top level

    white_mesa_network = whte_mesa.network = Network() # Network


    #reactor = SMPWR()  # Create reactor module
    #reactor.name = "SM-PWR"

    #reactor.shutdown = (True, 60*unit.minute)

    #white_mesa_network.module(reactor)  # Add reactor module to network

    # Leaching

    leaching = Leaching() # Create a Leaching module

    leaching.name = 'Leaching'
    leaching.save = True
    leaching.time_step = time_step
    leaching.end_time = end_time
    leaching.show_time = show_time

    white_mesa_network.module(leaching)

    # Decantation-Filtration

    decant_filt = DecantationFiltration()  # Create decantation/filtration module

    decant_filt.name = 'Decantation-Filtration'
    decant_filt.save = True
    decant_filt.time_step = time_step
    decant_filt.end_time = end_time
    decant_filt.show_time = show_time

    white_mesa_network.module(decant_filt)  # Add filtration module to network

    # Solvent extraction

    solvex = Solvex()  # Create solvent extraction module

    solvex.name = 'Solvex'
    solvex.save = True
    solvex.time_step = time_step
    solvex.end_time = end_time
    solvex.show_time = show_time

    white_mesa_network.module(solvex)  # Add solvex module to network

    # Precipitation

    precipt = Precipitation()  # Create solvent extraction module

    precipt.name = 'Precipitation'
    precipt.save = True
    precipt.time_step = time_step
    precipt.end_time = end_time
    precipt.show_time = show_time

    white_mesa_network.module(precipt)  # Add solvex module to network

    # Evaporation-Calcination

    evap_calc = EvaporationCalcination()  # Create solvent extraction module

    evap_calc.name = 'Evaporation-Calcination'
    evap_calc.save = True
    evap_calc.time_step = time_step
    evap_calc.end_time = end_time
    evap_calc.show_time = show_time

    white_mesa_network.module(evap_calc)  # Add solvex module to network

    # White Mesa Network Connectivity

    white_mesa_network.connect([leaching, 'pre-leach-product'], [decant_filt, 'std-feed'])
    white_mesa_network.connect([decant_filt, 'ccd-overflow'], [leaching, 'pre-leach-feed'])
    white_mesa_network.connect([leaching, 'acid-leach-product'], [decant_filt, 'ccd-feed'])
    white_mesa_network.connect([decant_filt, 'std-underflow'], [leaching, 'acid-leach-feed'])
    white_mesa_network.connect([decant_filt, 'filtrate'], [solvex, 'extraction-feed'])
    white_mesa_network.connect([solvex, 'raffinate'], [decant_filt, 'raffinate-feed'])
    white_mesa_network.connect([solvex, 'product'], [precipt, 'uts-feed'])
    white_mesa_network.connect([precipt, 'adu-product'], [evap_calc, 'adu-feed'])

    white_mesa_network.draw(engine='dot', node_shape='folder', size='1200,600', lr=True)
    #white_mesa_network.draw(engine='osage', node_shape='folder', size='1200,600', 
                            #graph_attr={'splines':'true', 'overlap':'scale', 'ranksep':'1.8'})
    #                        graph_attr={'splines':'true', 'ranksep':'3.0'})
    #white_mesa_network.draw()

    # Run
    if make_run:

        for module in white_mesa_network.modules:
            module.time_step = time_step
            module.end_time = end_time
            module.show_time = show_time

        whte_mesa.run()  # Run network dynamics simulation

    # Cortix run closure
    whte_mesa.close()  # Properly shutdow whte_mesa

    # Plots
    if make_plots and whte_mesa.use_multiprocessing or whte_mesa.rank == 0:

        # Reactor plots
        reactor = white_mesa_network.modules[0]

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
        plt.savefig('reactor-flowrate.png', dpi=300)

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

        # Steamer plots
        steamer = white_mesa_network.modules[1]

        (quant, time_unit) = steamer.primary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-primary-outflow-temp.png', dpi=300)

        (quant, time_unit) = steamer.secondary_inflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-secondary-inflow-temp.png', dpi=300)

        (quant, time_unit) = steamer.secondary_inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('steamer-secondary-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = steamer.secondary_outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('steamer-secondary-outflow-temp.png', dpi=300)

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
        plt.savefig('steamer-nusselt_s.png', dpi=300)

        # Turbine plots
        turbine = white_mesa_network.modules[2]

        (quant, time_unit) = turbine.state_phase.get_quantity_history('power')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('turbine-power.png', dpi=300)

        (quant, time_unit) = turbine.outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('turbine-outflow-temp.png', dpi=300)

        (quant, time_unit) = turbine.state_phase.get_quantity_history('rejected-heat')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('turbine-rejected-heat.png', dpi=300)

        # Condenser plots
        condenser = white_mesa_network.modules[3]

        (quant, time_unit) = condenser.inflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('condenser-inflow-temp.png', dpi=300)

        (quant, time_unit) = condenser.inflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('condenser-inflow-flowrate.png', dpi=300)

        (quant, time_unit) = condenser.outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('condenser-outflow-temp.png', dpi=300)

        # Water heater plots
        water_heater = white_mesa_network.modules[4]

        (quant, time_unit) = water_heater.outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('water_heater-outflow-temp.png', dpi=300)

        (quant, time_unit) = water_heater.outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+r' ['+quant.unit+']')
        plt.grid()
        plt.savefig('water_heater-flowrate.png', dpi=300)

        (quant, time_unit) = water_heater.inflow_phase.get_quantity_history('external-heat')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+r' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('water_heater-external-heat.png', dpi=300)

        (quant, time_unit) = water_heater.outflow_phase.get_quantity_history('rejected-heat')

        quant.plot(x_scaling=1/unit.minute, y_scaling=1/unit.mega, x_label='Time [m]',
                   y_label=quant.latex_name+r' [M'+quant.unit+']')
        plt.grid()
        plt.savefig('water_heater-rejected-heat.png', dpi=300)

if __name__ == '__main__':
    main()
