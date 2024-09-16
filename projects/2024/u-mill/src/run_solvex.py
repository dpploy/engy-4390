#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

from cortix import Cortix
from cortix import Network
from cortix import Units as unit

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

if __name__ == '__main__':
    main()
