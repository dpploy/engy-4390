#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""
Cortix Module
This module is a model of the Evaporation/Calcination process


             Ammonium Diuranate Feed
                       |
                       |
                       v
             ____________________
             |                  |
             |   Evaporation    |<----- Evaporator Heat
             |                  |
             |   Calcination    |<----- Calciner Heat
             |__________________|
                |            |
                |            |
                v            v
             Product      Off-Gas
              (U3O8)

 + Evaporation:
     1) Evaporator Base Parameters
       - # of Evaporator Columsns:
       - volume per Evaporator:                 1.42 m^3
       - Temperature Setpoint per Evaporator:   351.85 C
       - Feed Flowrate Entering Evaporator:    0.09 gallons/min
       - feed mass fraction of solids:
       - wash water mass fraction of solids:
       - overflow mass fraction of solids:
       - underflow mass fraction of solids:

 + Calcination:
     1) Calciner Heat Parameters
       - # of Rotary Calciners:                    2
       - volume per Calciner:                    m^3
       - Temperature Setpoint per Calciner:     850.0 C
       - Feed Flowrate Entering Evaporator:    0.07 gallons/min
       - feed mass fraction of solids:
       - wash water mass fraction of solids:
       - overflow mass fraction of solids:
       - underflow mass fraction of solids:

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

class EvaporationCalcination(Module):
    """Evaporation-calcination system.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: Precipitation, .
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['adu-feed', 'product']

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
        self.diuranate_volume = 0.988 * unit.meter**3
        self.entering_DUO_flowrate = 0.04 * unit.gallon/unit.minute
        self.calcination_flash_time = 80 * unit.minute

        #Evaporation

        #Calcination
        self.evaporation_feed_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.evaporation_feed_mass_density = 7.8 * unit.kg/unit.liter
        self.evaporation_feed_solids_massfrac = 100 * unit.ppm

        self.calcination_product_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.calcination_product_mass_density = 7.8 * unit.kg/unit.liter
        self.calcination_product_solids_massfrac = 100 * unit.ppm

        # Initialization
        '''
        self.primary_inflow_temp = primary_inflow_temp

        self.primary_pressure = 127.6*unit.bar

        self.primary_mass_flowrate = 4.66e6*unit.lb/unit.hour

        self.primary_outflow_temp = self.primary_inflow_temp #- 2*unit.K
        '''

        # Derived quantities
        '''
        self.rho_p = 0.0
        self.cp_p = 0.0
        self.mu_p = 0.0
        self.k_p = 0.0
        self.prtl_p = 0.0
        self.rey_p = 0.0
        self.nusselt_p = 0.0

        self.heat_sink_pwr = 0.0
        '''
        #***************************************************************************************
        # E V A P O R A T I O N
        #***************************************************************************************

        # Evaporation Feed Phase History (precipitation outflow)
        quantities = list()
        species = list()

        feed_mass_flowrate = Quantity(name='mass_flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.evaporation_feed_mass_flowrate,
                        latex_name=r'$\dot{m}$',
                        info='Evaporation Feed Mass Flowrate')
        quantities.append(feed_mass_flowrate)

        feed_mass_density = Quantity(name='mass_density',
                        formal_name='rho', unit='kg/m^3',
                        value=self.evaporation_feed_mass_density,
                        latex_name=r'$\rho$',
                        info='Evaporation Feed Mass Density')
        quantities.append(feed_mass_density)

        feed_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.evaporation_feed_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Evaporation Feed Solids Mass Fraction')

        quantities.append(feed_solids_massfrac)

        diuranate_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(s)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(diuranate_feed)

        self.evaporation_feed_phase = Phase(time_stamp=self.initial_time,
                                            time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # C A L C I N A T I O N
        #***************************************************************************************

        # Calcination Product Phase History (precipitation outflow)
        quantities = list()
        species = list()

        product_mass_flowrate = Quantity(name='mass_flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.calcination_product_mass_flowrate,
                        latex_name=r'$\dot{m}$',
                        info='Calcination Product Mass Flowrate')
        quantities.append(feed_mass_flowrate)

        product_mass_density = Quantity(name='mass_density',
                        formal_name='rho', unit='kg/m^3',
                        value=self.calcination_product_mass_density,
                        latex_name=r'$\rho$',
                        info='Calcination Product Mass Density')
        quantities.append(product_mass_density)

        product_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.calcination_product_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Calcination Product Solids Mass Fraction')

        quantities.append(product_solids_massfrac)

        u3o8_product = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(s)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(u3o8_product)

        self.calcination_product_phase = Phase(time_stamp=self.initial_time,
                                               time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # S T A T E  P H A S E
        #***************************************************************************************

        '''
        quantities = list()

        tau_p = Quantity(name='tau_p',
                        formal_name='Tau_p', unit='s',
                        value=0.0,
                        latex_name=r'$\tau_{p}$',
                        info='Steamer Primary Residence Time')

        quantities.append(tau_p)

        tau_s = Quantity(name='tau_s',
                        formal_name='Tau_s', unit='s',
                        value=0.0,
                        latex_name=r'$\tau_{s}$',
                        info='Steamer Secondary Residence Time')

        quantities.append(tau_s)

        heatflux = Quantity(name='heatflux',
                        formal_name="q''", unit='W/m$^2$',
                        value=0.0,
                        latex_name=r"$q''$",
                        info='Steamer Heat Flux')

        quantities.append(heatflux)

        nusselt_p = Quantity(name='nusselt_p',
                        formal_name='Nu_p', unit='',
                        value=0.0,
                        latex_name=r'$Nu_p$',
                        info='Steamer Primary Nusselt Number')

        quantities.append(nusselt_p)

        nusselt_s = Quantity(name='nusselt_s',
                        formal_name='Nu_s', unit='',
                        value=0.0,
                        latex_name=r'$Nu_s$',
                        info='Steamer Secondary Nusselt Number')

        quantities.append(nusselt_s)

        self.state_phase = Phase(time_stamp=self.initial_time,
                                 time_unit='s', quantities=quantities)
        '''

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

        # Interactions in the primary-inflow port
        #----------------------------------------
        # One way "from" adu-feed

        # Receive from
        if self.get_port('adu-feed').connected_port:

            self.send(time, 'adu-feed')

            (check_time, adu_feed) = self.recv('adu-feed')
            assert abs(check_time-time) <= 1e-6

            '''
            self.primary_inflow_temp = primary_inflow['temperature']
            self.primary_ressure = primary_inflow['pressure']
            self.primary_mass_flowrate = primary_inflow['mass_flowrate']
            '''

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # One way "to" product

        # Send to
        if self.get_port('product').connected_port:

            msg_time = self.recv('product')

            #temp = self.calcination_product_phase.get_value('temp', msg_time)

            product = dict()
            '''
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_pressure
            primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
            primary_outflow['quality'] = 0.0
            '''

            self.send((msg_time, product), 'product')

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
