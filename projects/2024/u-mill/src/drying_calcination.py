#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module.
   Drying/Calcination process in the White Mesa Uranium Milling Plant.


            Ammonium Diuranate Feed
                     |
                     |
                     v
             ____________________
 Off-Gas     |                  |
    <--------|     Drying       |<----- Steam Sparging (internal)
             |                  |<----- Resistance Heating (internal)
             |------------------|
             |                  |<----- Sweeping Gas
             |   Calcination    |<----- Resistance Heating (internal)
             |__________________|
                |            |
                |            |
                v            v
             Product      Off-Gas
              (U3O8)
+ Drying:
     1) Drying Base Parameters
       - 6 hearth furnace
       - # of Drying Columns:                 5
       - Volume per Drying:                 200 m^3
       - Temperature Setpoint per Drying:   351.85 C
       - Feed Flowrate Entering Drying:    0.09 gallons/min
       - feed mass fraction of solids:
       - wash water mass fraction of solids:
       - overflow mass fraction of solids:
       - underflow mass fraction of solids:

+ Calcination:
     1) Calciner Heat Parameters
       - # of rotary calciners:                    2
       - 6 hearth furnace
       - Volume per calciner:                    350 m^3
       - Temperature Setpoint per Calciner:     850.0 C
       - Feed flowrate entering calciner:    0.07 gallons/min
       - Feed mass fraction of solids:
       - Wash water mass fraction of solids:
       - Overflow mass fraction of solids:
       - Underflow mass fraction of solids:

Source of info:
 1993 Uranium Extraction Technology, IAEA Technical Reports Series No. 359
  p. ??? (White Mesa),
  p. ??? (Drying/Calc)
 URL: https://www-pub.iaea.org/MTCD/Publications/PDF/trs359_web.pdf

"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity
from cortix import Species
from cortix import Units as unit

class DryingCalcination(Module):
    """Drying-calcination system.

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

        # Drying
        self.dryer_tank_volume = 2 * 200 * unit.meter**3
        self.five_tau_drying = 48 * unit.hour
        self.dry_densification_factor = 1.2
        self.dry_shrink_factor = 1.1
        self.dry_volume_reduction = 50/100
        self.product_discharge_time = 1.0 * unit.day
        self.dryer_status = 'empty' # filling-up, drying, discharging

        # Calcination

        # Initialization

        # Drying
        if self.get_port('adu-feed').connected_port:
            self.dry_feed_mass_flowrate = 0 * unit.kg/unit.minute
            self.dry_feed_mass_density = 0 * unit.kg/unit.liter
            self.dry_feed_solids_massfrac = 0 * unit.ppm
        else:
            self.dry_feed_mass_flowrate = 5800.0 * unit.kg/unit.minute
            self.dry_feed_mass_density = 8.7 * unit.kg/unit.liter
            self.dry_feed_solids_massfrac = 100 * unit.ppm

        # Calcination

        # Derived quantities

        #***************************************************************************************
        # E V A P O R A T I O N
        #***************************************************************************************

        # Drying Feed Phase History (precipitation outflow)
        quantities = list()
        species = list()

        mass_flowrate = Quantity(name='mass-flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=0.0,
                        latex_name=r'$\dot{m}$',
                        info='Drying Product Mass Flowrate')
        quantities.append(mass_flowrate)

        mass_density = Quantity(name='mass-density',
                        formal_name='rho', unit='kg/m^3',
                        value=0.0,
                        latex_name=r'$\rho$',
                        info='Drying Product Mass Density')
        quantities.append(mass_density)

        solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=0.0,
                        latex_name=r'$C_1$',
                        info='Drying Product Solids Mass Fraction')

        quantities.append(solids_massfrac)

        diuranate = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(s)',
                            atoms=['U','2*O','3*S','12*O'],
                            info='UO2-(SO4)3^4-')
        species.append(diuranate)

        self.dry_product_phase = Phase(time_stamp=self.initial_time,
                                        time_unit='s', quantities=quantities, species=species)

        # Drying Off-Gas Phase History
        quantities = list()
        species = list()

        mass_flowrate = Quantity(name='mass-flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=0.0,
                        latex_name=r'$\dot{m}$',
                        info='Drying Off-Gas Mass Flowrate')
        quantities.append(mass_flowrate)

        mass_density = Quantity(name='mass-density',
                        formal_name='rho', unit='kg/m^3',
                        value=0.0,
                        latex_name=r'$\rho$',
                        info='Drying Off-Gas Mass Density')
        quantities.append(mass_density)

        solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=0.0,
                        latex_name=r'$C_1$',
                        info='Drying Off-Gas Solids Mass Fraction')

        quantities.append(solids_massfrac)

        diuranate = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(s)',
                            atoms=['U','2*O','3*S','12*O'],
                            info='UO2-(SO4)3^4-')
        species.append(diuranate)

        self.dry_offgas_phase = Phase(time_stamp=self.initial_time,
                                       time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # C A L C I N A T I O N
        #***************************************************************************************

        # Calcination Product Phase History (precipitation outflow)
        quantities = list()
        species = list()

        product_mass_flowrate = Quantity(name='mass-flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=0.0,
                        latex_name=r'$\dot{m}$',
                        info='Calcination Product Mass Flowrate')
        quantities.append(product_mass_flowrate)

        product_mass_density = Quantity(name='mass_density',
                        formal_name='rho', unit='kg/m^3',
                        value=0.0,
                        latex_name=r'$\rho$',
                        info='Calcination Product Mass Density')
        quantities.append(product_mass_density)

        product_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=0.0,
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

        # Drying
        quantities = list()
        species = list()

        liq_volume = Quantity(name='liquid-volume',
                        formal_name='v', unit='m$^3$',
                        value=0.0,
                        latex_name=r'$V_{e}$',
                        info='Drying Tank Liquid Volume')
        quantities.append(liq_volume)

        self.dry_state_phase = Phase(time_stamp=self.initial_time, time_unit='s',
                                      quantities=quantities, species=species)

    def run(self, *args):

        # Some logic for logging time stamps
        # Leave this here: rebuild logger
        logger_name = args[0][0].name
        self.rebuild_logger(logger_name)

        self.end_time = max(self.end_time, self.initial_time + self.time_step)

        time = self.initial_time

        print_time = self.initial_time
        print_time_step = self.show_time[1]

        print_time_step = max(print_time_step, self.time_step)

        while time <= self.end_time:

            if self.show_time[0] and \
               (print_time <= time < print_time+print_time_step):

                msg = self.name+'::run():time[d]='+ str(round(time/unit.day, 1))
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

        #-----------------------------
        # Evolve the Drying state
        #-----------------------------
        rho_dry_ideal = self.dry_feed_mass_density * self.dry_densification_factor

        if self.dryer_status == 'empty' or self.dryer_status == 'filling-up':
            rho_dry_feed = self.dry_feed_mass_density
            mass_flowrate_dry_feed = self.dry_feed_mass_flowrate
            vol_flowrate_feed = mass_flowrate_dry_feed/rho_dry_feed

            rho_dry = rho_dry_feed
            mass_flowrate_dry = 0
        else:
            rho_dry_feed = 0.0
            mass_flowrate_dry_feed = 0.0
            vol_flowrate_feed = 0.0

        dry_liq_volume_initial = self.dry_state_phase.get_value('liquid-volume', time)

        rho_dry_initial = self.dry_product_phase.get_value('mass-density', time)

        if self.dryer_status != 'discharging':
            vol_flowrate_dry_initial = 0.0
            mass_flowrate_dry = 0
        else:
            vol_flowrate_dry_initial = dry_liq_volume_initial/self.product_discharge_time
            mass_flowrate_dry = rho_dry_initial * vol_flowrate_dry_initial
            rho_dry = dry_liq_volume_initial

        dry_liq_volume = dry_liq_volume_initial + \
                          (vol_flowrate_feed - vol_flowrate_dry_initial) * self.time_step

        assert dry_liq_volume >= 0.0

        # Ideal Drying

        # Filling-up
        if dry_liq_volume < self.dryer_tank_volume and self.dryer_status == 'empty':
            self.dryer_status = 'filling-up'

        # Drying
        if (dry_liq_volume > self.dryer_tank_volume and self.dryer_status == 'filling-up') or \
           self.dryer_status == 'drying':
            self.dryer_status = 'drying'
            tau = self.five_tau_drying
            rho_dry = rho_dry_ideal + math.exp(-self.time_step/tau) * (rho_dry_initial - rho_dry_ideal)
            mass_flowrate_dry = 0.0

            delta_rho = rho_dry - rho_dry_initial
            delta_vol = self.dry_shrink_factor / delta_rho
            dry_liq_volume = dry_liq_volume_initial - delta_vol

            if dry_liq_volume <= self.dry_volume_reduction * self.dryer_tank_volume:
                self.dryer_status = 'discharging'
                print('Discharging')

        # Discharge
        if rho_dry_initial >= 0.95*rho_dry_ideal:
            self.dryer_status = 'discharging'

        tmp_dry_state = self.dry_state_phase.get_row(time)
        tmp_dry_product = self.dry_product_phase.get_row(time)


        '''mass_flowrate_initial = self.drying_phase.get_value('mass_flowrate', time)

        mass_flowrate_inflow = self.drying_feed_mass_flowrate

        rho_drying = self.drying_feed_mass_density

        vol_flowrate_initial = mass_flowrate_initial/rho_drying

        if vol_flowrate_initial == 0:
            vol_flowrate_initial = mass_flowrate_inflow
            tau = self.dryer_volume/vol_flowrate_initial
        else:
            tau = self.dryer_volume/vol_flowrate_initial

        # Mass balance
        mass_flowrate_drying = mass_flowrate_inflow + \
                                 math.exp(-time/tau) * (mass_flowrate_initial - mass_flowrate_inflow)
'''
        '''
        # Evolve the u3o8 Product state
        mass_flowrate_initial = self.calcination_product_phase.get_value('mass_flowrate', time)

        mass_flowrate_outflow = self.calcination_product_mass_flowrate

        rho_calcination = self.calcination_product_mass_density

        vol_flowrate_initial = mass_flowrate_initial/rho_calcination

        if vol_flowrate_initial == 0:
            vol_flowrate_initial = mass_flowrate_outflow
            tau = self.dryer_volume/vol_flowrate_initial
        else:
            tau = self.dryer_volume/vol_flowrate_initial

        # Mass balance
        mass_flowrate_product = mass_flowrate_outflow + \
                                  math.exp(-time/tau) * (mass_flowrate_initial - mass_flowrate_outflow)

        tmp_calcination = self.calcination_phase.get_row(time)
        '''

        #----------------------------
        # Step All Quantities in Time
        #----------------------------

        time += self.time_step

        # Drying
        self.dry_state_phase.add_row(time, tmp_dry_state)
        self.dry_state_phase.set_value('liquid-volume', dry_liq_volume, time)

        self.dry_product_phase.add_row(time, tmp_dry_product)
        self.dry_product_phase.set_value('mass-flowrate', mass_flowrate_dry, time)
        self.dry_product_phase.set_value('mass-density', rho_dry, time)

        '''
        # Calcination
        self.calcination_product_phase.add_row(time, tmp_calcination)
        #self.drying_phase.add_row(time, tmp_calcination)

        self.calcination_product_phase.set_value('mass_flowrate', mass_flowrate_product, time)
        self.calcination_product_phase.set_value('mass_density', rho_calcination, time)

        #self.drying_phase.set_value('mass_flowrate', mass_flowrate_drying, time)
        '''

        return time
