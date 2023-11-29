#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module.
   Decantation/Filtration process in the White Mesa Uranium Milling Plant.
   Add info here:...

   Mass fraction in ppm.

   + Decantation:
       - # of thickeners:                    4
       - volume per thickener:               0.988 m^3
       - wash water volumetric flowrate:     1067 gallons/min
       - feed mass fraction of solids:       100000 ppm
       - wash water mass fraction of solids: 500 ppm
       - overflow mass fraction of solids:   100 ppm
       - underflow mass fraction of solids:  99000 ppm

   + Filtration:

   Source of info:
"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity
from cortix import Species

import unit

class DecantationFiltration(Module):
    """Decantation-filtration system.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: Leaching, Solvex
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['feed', 'filtrate']

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)
        self.save = True

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters
        self.thickener_volume = 0.988 * unit.meter**3
        self.wash_water_vol_flowrate = 1067 * unit.gallon/unit.minute
        self.decantation_relaxation_time = 20 * unit.minute

        # Initialization

        # Decantation
        self.decantation_feed_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.decantation_feed_mass_density = 7.8 * unit.kg/unit.liter
        self.decantation_feed_solids_massfrac = 100 * unit.ppm

        self.decantation_underflow_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.decantation_underflow_solids_massfrac = 100 * unit.ppm

        self.decantation_overflow_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.decantation_overflow_solids_massfrac = 100 * unit.ppm

        # Filtration
        self.filtration_slurry_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_slurry_solids_massfrac = 100 * unit.ppm

        self.filtration_filtrate_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_filtrate_solids_massfrac = 100 * unit.ppm

        # Derived quantities
        # self.xxxx

        #***************************************************************************************
        # D E C A N T A T I O N
        #***************************************************************************************

        # Decantation Feed Phase History (leaching outflow)
        quantities = list()
        species = list()

        feed_mass_flowrate = Quantity(name='mass_flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.decantation_feed_mass_flowrate,
                        latex_name=r'$\dot{m}$',
                        info='Decantation Feed Mass Flowrate')
        quantities.append(feed_mass_flowrate)

        feed_mass_density = Quantity(name='mass_density',
                        formal_name='rho', unit='kg/m^3',
                        value=self.decantation_feed_mass_density,
                        latex_name=r'$\rho$',
                        info='Decantation Feed Mass Density')
        quantities.append(feed_mass_density)

        feed_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.decantation_feed_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Decantation Feed Solids Mass Fraction')

        quantities.append(feed_solids_massfrac)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(aq)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O',formula_name='H2O(aq)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        h2so4_feed = Species(name='H2SO4',formula_name='H2SO4(aq)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_feed)

        iron_feed = Species(name='Fe', formula_name='Fe(s)',
                           atoms=['Fe'],
                           info='Fe')
        species.append(iron_feed)

        copper_feed = Species(name='Cu', formula_name='Cu(s)',
                           atoms=['Cu'],
                           info='Cu')
        species.append(copper_feed)

        gold_feed = Species(name='Au', formula_name='Au(s)',
                           atoms=['Au'],
                           info='Au')
        species.append(gold_feed)

        self.decantation_feed_phase = Phase(time_stamp=self.initial_time,
                                            time_unit='s', quantities=quantities, species=species)

        # Decantation Underflow Phase History (internal state/external)
        quantities = list()
        species = list()  # Is it proper to rest the species list too?

        underflow_mass_flowrate = Quantity(name='mass_flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=self.decantation_underflow_mass_flowrate,
                                       latex_name=r'$\dot{m}_4$',
                                       info='Decantation Underflow Mass Flowrate')
        quantities.append(underflow_mass_flowrate)

        underflow_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.decantation_underflow_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Decantation Underflow Solids Mass Fraction')

        quantities.append(underflow_solids_massfrac)

        uo2so434minus_underflow = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(aq)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_underflow)

        h2o_underflow = Species(name='H2O',formula_name='H2O(aq)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_underflow)

        h2so4_underflow = Species(name='H2SO4',formula_name='H2SO4(aq)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_underflow)

        iron_underflow = Species(name='Fe', formula_name='Fe(s)',
                           atoms=['Fe'],
                           info='Fe')
        species.append(iron_underflow)

        copper_underflow = Species(name='Cu', formula_name='Cu(s)',
                           atoms=['Cu'],
                           info='Cu')
        species.append(copper_underflow)

        gold_underflow = Species(name='Au', formula_name='Au(s)',
                           atoms=['Au'],
                           info='Au')
        species.append(gold_underflow)

        self.decantation_underflow_phase = Phase(time_stamp=self.initial_time,
                                                 time_unit='s', quantities=quantities, species=species)

        # Decantation Overflow Phase History (internal state)
        quantities = list()
        species = list()

        overflow_mass_flowrate = Quantity(name='mass_flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.decantation_overflow_mass_flowrate,
                        latex_name=r'$\dot{m}_1$',
                        info='Decantation Overflow Mass Flowrate')
        quantities.append(overflow_mass_flowrate)

        overflow_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.decantation_overflow_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Decantation Overflow Solids Mass Fraction')

        quantities.append(overflow_solids_massfrac)

        uo2so434minus_overflow = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(aq)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_overflow)

        h2o_overflow = Species(name='H2O',formula_name='H2O(aq)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_overflow)

        h2so4_overflow = Species(name='H2SO4',formula_name='H2SO4(aq)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_overflow)

        iron_overflow = Species(name='Fe', formula_name='Fe(s)',
                           atoms=['Fe'],
                           info='Fe')
        species.append(iron_overflow)

        copper_overflow = Species(name='Cu', formula_name='Cu(s)',
                           atoms=['Cu'],
                           info='Cu')
        species.append(copper_overflow)

        gold_overflow = Species(name='Au', formula_name='Au(s)',
                           atoms=['Au'],
                           info='Au')
        species.append(gold_overflow)

        self.decantation_overflow_phase = Phase(time_stamp=self.initial_time,
                                                time_unit='s', quantities=quantities, species=species)


        #***************************************************************************************
        # F I L T R A T I O N
        #***************************************************************************************

        # Filtration Slurry Phase History (internal state/external)
        quantities = list()
        species = list()  # Is it proper to rest the species list too?

        slurry_mass_flowrate = Quantity(name='mass_flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=self.filtration_slurry_mass_flowrate,
                                       latex_name=r'$\dot{m}_4$',
                                       info='Filtration Slurry Mass Flowrate')
        quantities.append(slurry_mass_flowrate)

        slurry_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.filtration_slurry_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Filtration Slurry Solids Mass Fraction')

        quantities.append(slurry_solids_massfrac)

        uo2so434minus_slurry = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(aq)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_slurry)

        h2o_slurry = Species(name='H2O',formula_name='H2O(aq)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_slurry)

        h2so4_slurry = Species(name='H2SO4',formula_name='H2SO4(aq)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_slurry)

        iron_slurry = Species(name='Fe', formula_name='Fe(s)',
                           atoms=['Fe'],
                           info='Fe')
        species.append(iron_slurry)

        copper_slurry = Species(name='Cu', formula_name='Cu(s)',
                           atoms=['Cu'],
                           info='Cu')
        species.append(copper_slurry)

        gold_slurry = Species(name='Au', formula_name='Au(s)',
                           atoms=['Au'],
                           info='Au')
        species.append(gold_slurry)

        self.filtration_slurry_phase = Phase(time_stamp=self.initial_time,
                                  time_unit='s', quantities=quantities, species=species)

        # Filtration Filtrate Phase History (internal state/external)
        quantities = list()
        species = list()  # Is it proper to rest the species list too?

        filtrate_mass_flowrate = Quantity(name='mass_flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=self.filtration_filtrate_mass_flowrate,
                                       latex_name=r'$\dot{m}_4$',
                                       info='Filtration Filtrate Mass Flowrate')
        quantities.append(filtrate_mass_flowrate)

        filtrate_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.filtration_filtrate_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Filtration Filtrate Solids Mass Fraction')

        quantities.append(filtrate_solids_massfrac)

        uo2so434minus_filtrate = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(aq)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_filtrate)

        h2o_filtrate = Species(name='H2O',formula_name='H2O(aq)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_filtrate)

        h2so4_filtrate = Species(name='H2SO4',formula_name='H2SO4(aq)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_filtrate)

        iron_filtrate = Species(name='Fe', formula_name='Fe(s)',
                           atoms=['Fe'],
                           info='Fe')
        species.append(iron_filtrate)

        copper_filtrate = Species(name='Cu', formula_name='Cu(s)',
                           atoms=['Cu'],
                           info='Cu')
        species.append(copper_filtrate)

        gold_filtrate = Species(name='Au', formula_name='Au(s)',
                           atoms=['Au'],
                           info='Au')
        species.append(gold_filtrate)

        self.filtration_filtrate_phase = Phase(time_stamp=self.initial_time,
                                  time_unit='s', quantities=quantities, species=species)

    def run(self, *args):

        # Some logic for logging time stamps
        if self.initial_time + self.time_step > self.end_time:
            self.end_time = self.initial_time + self.time_step

        time = self.initial_time

        print_time = self.initial_time
        print_time_step = self.show_time[1]

        if print_time_step < self.time_step:
            print_time_step = self.time_step

        while time <= self.end_time:

            if self.show_time[0] and \
               (print_time <= time < print_time+print_time_step):

                msg = self.name+'::run():time[m]='+ str(round(time/unit.minute, 1))
                self.log.info(msg)

                self.__logit = True
                print_time += self.show_time[1]

            else:
                self.__logit = False

            # Evolve one time step
            #---------------------

            time = self.__step(time)

            # Communicate information
            #------------------------

            self.__call_ports(time)

        self.end_time = time # correct the final time if needed

    def __call_ports(self, time):

        # Interactions in the feed port
        #------------------------------
        # One way "from" feed port

        # Receive from
        if self.get_port('feed').connected_port:

            self.send(time, 'feed')

            (check_time, feed_phase) = self.recv('feed')
            assert abs(check_time-time) <= 1e-6

            # insert data from feed_phase into decantation_feed_phase history

        # Interactions in the filtrate port
        #----------------------------------
        # One way "to" filtrate port

        # Send to
        if self.get_port('filtrate').connected_port:

            msg_time = self.recv('filtrate')

            '''
            temp = self.primary_outflow_phase.get_value('temp', msg_time)
            primary_outflow = dict()
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_pressure
            primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
            primary_outflow['quality'] = 0.0
            '''
            # extract filtration_filtrate_phase data at msg_time to send 

            self.send((msg_time, filtrate_phase), 'filtrate')

    def __step(self, time=0.0):
        """Stepping Decantation-Filtration in time
        """

        '''
        # Get state values
        u_0 = self.__get_state_vector(time)

        t_interval_sec = np.linspace(time, time+self.time_step, num=2)

        max_n_steps_per_time_step = 1500 # max number of nonlinear algebraic solver
                                         # iterations per time step

        (u_vec_hist, info_dict) = odeint(self.__f_vec, u_0, t_interval_sec,
                                         rtol=1e-7, atol=1e-8,
                                         mxstep=max_n_steps_per_time_step,
                                         full_output=True, tfirst=False)

        assert info_dict['message'] == 'Integration successful.', info_dict['message']

        u_vec = u_vec_hist[1, :]  # solution vector at final time step

        temp_p = u_vec[0] # primary outflow temp
        temp_s = u_vec[1] # secondary outflow temp

        # Update phases
        primary_outflow = self.primary_outflow_phase.get_row(time)
        secondary_inflow = self.secondary_inflow_phase.get_row(time)
        secondary_outflow = self.secondary_outflow_phase.get_row(time)
        steamer = self.state_phase.get_row(time)
        '''

        time += self.time_step

        '''
        self.primary_outflow_phase.add_row(time, primary_outflow)
        self.primary_outflow_phase.set_value('temp', temp_p, time)
        self.primary_outflow_phase.set_value('flowrate', self.primary_mass_flowrate, time)

        self.secondary_inflow_phase.add_row(time, secondary_inflow)
        self.secondary_inflow_phase.set_value('temp', self.secondary_inflow_temp, time)
        self.secondary_inflow_phase.set_value('flowrate', self.secondary_mass_flowrate, time)

        self.secondary_outflow_phase.add_row(time, secondary_outflow)
        self.secondary_outflow_phase.set_value('temp', temp_s, time)
        self.secondary_outflow_phase.set_value('flowrate', self.secondary_mass_flowrate, time)
        self.secondary_outflow_phase.set_value('pressure', self.secondary_pressure, time)
        self.secondary_outflow_phase.set_value('quality', self.secondary_outflow_quality, time)

        self.state_phase.add_row(time, steamer)

        # Primary residence time
        self.state_phase.set_value('tau_p', self.tau_p, time)

        # Secondary residence time
        self.state_phase.set_value('tau_s', self.tau_s, time)

        # Heat flux and Nusselt number
        heatflux = -self.heat_sink_pwr/self.heat_transfer_area
        self.state_phase.set_value('heatflux', heatflux, time)

        self.state_phase.set_value('nusselt_p', self.nusselt_p, time)

        self.state_phase.set_value('nusselt_s', self.nusselt_s, time)
        '''

        return time

