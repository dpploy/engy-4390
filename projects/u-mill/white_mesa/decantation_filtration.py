#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module.
   Decantation/Filtration process in the White Mesa Uranium Milling Plant.


             |-----------------|
             |                 |
             |   DECANTATION   |
 STD         |                 |
    <--------|  + Single Tank  |<-------- Wash water (internal)
 underflow   |                 |<-------- Pre-Leach Feed
             |                 |
 CCD         |                 |<-------- Acid-Leach Feed
    <--------|  + CC Bank Tank |<-------- Wash water (internal)
 overflow    |                 |<-------- Raffinate Feed (from solvent extraction)
             |-----------------|
             |                 |
             |   FILTRATION    |<-------- STD Overflow (internal)
             |                 |
             |-----------------|
                  |       |
                  |       |
                  |       |
                  v       v
               Filtrate  Slurry Waste (TBD)

 NB. STD underflow goes to acid leaching
 NB. CCD overflow goes to pre-leaching
 NB. Filtrate goes to solvent extraction


   + Decantation Steady-State Operation:
     1) Single-Tank Decantation  (STD)
       - volume of thickener:                3402.33 m^3
       - preleach feed mass flowrate:        4540 kg/hr
       - wash water volumetric flowrate:     118735 gallons/min or 99(feed flow rate)
       - overflow volumetric flowrate:       38 gallons/min
       - overflow mass fraction of solids:   100 ppm
       - feed mass fraction of solids:       10000 ppm
       - wash water mass fraction of solids: 0 ppm
       - underflow mass fraction of solids:  9900 ppm
       - wash water ratio (R):               99

     2) Counter-Current Decantation (CCD)
       - # of thickeners:                    7
       - volume per thickener:               339.29 m^3
       - acid leach feed mass flowrate:
       - wash water volumetric flowrate:     gallons/min or 1.69(feed flow rate)
       - feed mass fraction of solids:       10000 ppm
       - wash water mass fraction of solids: 500 ppm
       - overflow mass fraction of solids:   100 ppm
       - underflow mass fraction of solids:  9900 ppm
       - wash water ratio (R):               1.69

   + Filtration Steady-State Operation:
       - volume of Drum Filter:              32.04 m^3
       - filtrate mass fraction of solids:   10 ppm
       - slurry mass fraction of solids:     990 ppm

   Source of info:
   -https://www.911metallurgist.com/blog/uranium-extraction-process
   Filtration/Thickeners Section

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

        self.port_names_expected = ['std-feed', 'ccd-feed', 'raffinate-feed',
                                    'filtrate', 'std-underflow', 'ccd-overflow']

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
        self.std_tank_volume = 3402.33 * unit.meter**3
        self.wash_water_preleach_feed_ratio = 99.0
        self.std_underflow_overflow_mass_flowrate_ratio = 99.0
        self.wash_water_acidleach_feed_ratio = 1.69

        self.ccd_tank_volume = 7*339.29 * unit.meter**3 # CCD
        self.ccd_wash_water_vol_flowrate = 1067 * unit.gallon/unit.minute
        self.ccd_underflow_overflow_mass_flowrate_ratio = 99.0

        # Initialization
        self.wash_water_mass_density = 1.0 * unit.kg/unit.meter**3

        # Decantation
        if self.get_port('std-feed').connected_port:
            self.single_tank_decantation_preleach_feed_mass_flowrate = 0 * unit.kg/unit.minute
            self.single_tank_decantation_preleach_feed_mass_density = 0 * unit.kg/unit.liter
            self.single_tank_decantation_preleach_feed_solids_massfrac = 0 * unit.ppm
        else:
            self.single_tank_decantation_preleach_feed_mass_flowrate = 4540 * unit.kg/unit.minute
            self.single_tank_decantation_preleach_feed_mass_density = 7.8 * unit.kg/unit.liter
            self.single_tank_decantation_preleach_feed_solids_massfrac = 100 * unit.ppm

        self.ccd_acidleach_feed_mass_flowrate = 2270 * unit.kg/unit.minute
        self.ccd_acidleach_feed_mass_density = 7.8 * unit.kg/unit.liter
        self.ccd_acidleach_feed_solids_massfrac = 100 * unit.ppm

        self.single_tank_decantation_raffinate_feed_mass_flowrate = 1.0 * unit.kg/unit.minute
        self.single_tank_decantation_raffinate_feed_mass_density = 7.8 * unit.kg/unit.liter
        self.single_tank_decantation_raffinate_feed_solids_massfrac = 150 * unit.ppm

        self.ccd_underflow_mass_flowrate = 1.0 * unit.kg/unit.minute
        self.ccd_underflow_solids_massfrac = 100 * unit.ppm

        self.ccd_overflow_mass_flowrate = 1.0 * unit.kg/unit.minute
        self.ccd_overflow_solids_massfrac = 100 * unit.ppm

        # Filtration
        self.filtration_slurry_mass_flowrate = 1.0 * unit.kg/unit.minute
        self.filtration_slurry_solids_massfrac = 100 * unit.ppm

        self.filtration_filtrate_mass_flowrate = 1.0 * unit.kg/unit.minute
        self.filtration_filtrate_solids_massfrac = 10 * unit.ppm

        # Derived quantities
        # self.xxxx

        #***************************************************************************************
        # D E C A N T A T I O N
        #***************************************************************************************

        # Single Tank Decantation Overflow Phase History (Interal feed to Filtration)
        quantities = list()
        species = list()

        feed_mass_flowrate = Quantity(name='mass-flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=0.0,
                        latex_name=r'$\dot{m}$',
                        info='Decantation Single-Tank Overflow Mass Flowrate')
        quantities.append(feed_mass_flowrate)

        feed_mass_density = Quantity(name='mass_density',
                        formal_name='rho', unit='kg/m^3',
                        value=0.0,
                        latex_name=r'$\rho$',
                        info='Decantation Single-Tank Overflow Feed Mass Density')
        quantities.append(feed_mass_density)

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

        self.single_tank_decantation_overflow_phase = Phase(time_stamp=self.initial_time,
                                            time_unit='s', quantities=quantities, species=species)

        # Single Tank Decantation Underflow Phase History (Send to Acid-leaching)
        quantities = list()
        species = list()

        underflow_mass_flowrate = Quantity(name='mass-flowrate',
                                      formal_name='mdot', unit='kg/s',
                                      value=0.0,
                                      latex_name=r'$\dot{m}$',
                                      info='Decantation Single-Tank Underflow Mass Flowrate')
        quantities.append(underflow_mass_flowrate)

        underflow_solids_massfrac = Quantity(name='solids_massfrac',
                                             formal_name='solids_massfrac', unit='ppm',
                                             value=0.0,
                                             latex_name=r'$C_1$',
                                             info='STD Underflow Solids Mass Fraction')
        quantities.append(underflow_solids_massfrac)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(aq)',
                                     atoms=['U', '2*O', '3*S', '12*O'],
                                     info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O', formula_name='H2O(aq)',
                           atoms=['2*H', 'O'],
                           info='H2O')
        species.append(h2o_feed)

        h2so4_feed = Species(name='H2SO4', formula_name='H2SO4(aq)',
                             atoms=['2*H', 'S', '4*O'],
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

        self.single_tank_decantation_underflow_phase = Phase(time_stamp=self.initial_time,
                                                                  time_unit='s', quantities=quantities, species=species)

        # Counter Current Decantation Underflow Phase History (Tailings)
        quantities = list()
        species = list()

        ccd_underflow_mass_flowrate = Quantity(name='mass-flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=self.ccd_underflow_mass_flowrate,
                                       latex_name=r'$\dot{m}_4$',
                                       info='Counter-Current Decantation Underflow Mass Flowrate')
        quantities.append(ccd_underflow_mass_flowrate)

        underflow_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.ccd_underflow_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Counter-Current Decantation Underflow Solids Mass Fraction')

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

        self.ccd_underflow_phase = Phase(time_stamp=self.initial_time,
                                                 time_unit='s', quantities=quantities, species=species)

        # Counter Current Decantation Overflow Phase History (Send to Pre-leaching)
        quantities = list()
        species = list()

        overflow_mass_flowrate = Quantity(name='mass-flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.ccd_overflow_mass_flowrate,
                        latex_name=r'$\dot{m}_1$',
                        info='Counter-Current Decantation Overflow Mass Flowrate')
        quantities.append(overflow_mass_flowrate)

        overflow_solids_massfrac = Quantity(name='solids_massfrac',
                        formal_name='solids_massfrac', unit='ppm',
                        value=self.ccd_overflow_solids_massfrac,
                        latex_name=r'$C_1$',
                        info='Counter-Current Decantation Overflow Solids Mass Fraction')

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

        self.ccd_overflow_phase = Phase(time_stamp=self.initial_time,
                                                time_unit='s', quantities=quantities, species=species)


        #***************************************************************************************
        # F I L T R A T I O N
        #***************************************************************************************

        # Filtration Slurry Waste Phase History (Tailings/tbd)
        quantities = list()
        species = list()

        slurry_mass_flowrate = Quantity(name='mass-flowrate',
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

        # Filtration Filtrate Phase History (Send to Solvent Extraction)
        quantities = list()
        species = list()

        filtrate_mass_flowrate = Quantity(name='mass-flowrate',
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

        #***************************************************************************************
        # S T A T E  P H A S E
        #***************************************************************************************

        # Single Tank Decantation State
        quantities = list()
        species = list()

        liq_volume = Quantity(name='liquid-volume',
                        formal_name='v', unit='m$^3$',
                        value=0.0,
                        latex_name=r'$V_\text{std}$',
                        info='Decantation Single-Tank Liquid Volume')
        quantities.append(liq_volume)

        self.single_tank_decantation_state_phase = Phase(time_stamp=self.initial_time, time_unit='s',
                                                         quantities=quantities, species=species)

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

        # Interactions in the std feed port
        #------------------------------
        # One way "from" std feed port

        # Receive from
        if self.get_port('std-feed').connected_port:

            self.send(time, 'std-feed')

            (check_time, feed_stream) = self.recv('std-feed')
            assert abs(check_time-time) <= 1e-6

            # Update feed stream data
            self.single_tank_decantation_preleach_feed_mass_flowrate = feed_stream['mass-flowrate']
            self.single_tank_decantation_preleach_feed_mass_density = feed_stream['mass-density']
            self.single_tank_decantation_preleach_feed_solids_massfrac = feed_stream['solids-massfrac']

        # Interactions in the ccd feed port
        #------------------------------
        # One way "from" feed port

        # Receive from
        if self.get_port('ccd-feed').connected_port:

            self.send(time, 'ccd-feed')

            (check_time, feed_phase) = self.recv('ccd-feed')
            assert abs(check_time-time) <= 1e-6

            # insert data from ccd-feed_phase into decantation_feed_phase history

        # Interactions in the raffinate feed port
        #------------------------------
        # One way "from" raffinate feed port

        # Receive from
        if self.get_port('raffinate-feed').connected_port:

            self.send(time, 'raffinate-feed')

            (check_time, feed_phase) = self.recv('raffinate-feed')
            assert abs(check_time-time) <= 1e-6

            # insert data from raffinate-feed_phase into decantation_feed_phase history

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

            self.send((msg_time, self.filtration_filtrate_phase), 'filtrate')

            # Interactions in the ccd overflow port
            # ----------------------------------
            # One way "to" ccd overflow port

            # Send to
            if self.get_port('ccd-overflow').connected_port:
                msg_time = self.recv('ccd-overflow')

                '''
                temp = self.primary_outflow_phase.get_value('temp', msg_time)
                primary_outflow = dict()
                primary_outflow['temperature'] = temp
                primary_outflow['pressure'] = self.primary_pressure
                primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
                primary_outflow['quality'] = 0.0
                '''
                # extract filtration_filtrate_phase data at msg_time to send

                #this is not filtration filtrate phase v
                self.send((msg_time, self.filtration_filtrate_phase), 'ccd-overflow')

            # Interactions in the std-underflow port
            # ----------------------------------
            # One way "to" std underflow port

            # Send to
            if self.get_port('std-underflow').connected_port:
                msg_time = self.recv('std-underflow')

                '''
                temp = self.primary_outflow_phase.get_value('temp', msg_time)
                primary_outflow = dict()
                primary_outflow['temperature'] = temp
                primary_outflow['pressure'] = self.primary_pressure
                primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
                primary_outflow['quality'] = 0.0
                '''
                # extract filtration_filtrate_phase data at msg_time to send

                self.send((msg_time, self.filtration_filtrate_phase), 'std-underflow')

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

        #----------------------
        # Evolve the STD state
        #----------------------
        # Ideal flow mixing
        std_underflow_mass_flowrate_initial = self.single_tank_decantation_underflow_phase.get_value('mass-flowrate', time)
        std_overflow_mass_flowrate_initial = self.single_tank_decantation_overflow_phase.get_value('mass-flowrate', time)

        preleach_feed_mass_flowrate = self.single_tank_decantation_preleach_feed_mass_flowrate
        rho_preleach_feed = self.single_tank_decantation_preleach_feed_mass_density

        if preleach_feed_mass_flowrate == 0.0:
            preleach_feed_vol_flowrate = 0.0
        else:
            preleach_feed_vol_flowrate = preleach_feed_mass_flowrate/rho_preleach_feed

        wash_water_vol_flowrate = self.wash_water_preleach_feed_ratio * preleach_feed_vol_flowrate
        rho_wash_water = self.wash_water_mass_density
        wash_water_mass_flowrate = wash_water_vol_flowrate * rho_wash_water

        # Ideal solution
        mass_flowrate_inflow = preleach_feed_mass_flowrate + wash_water_mass_flowrate
        rho_std = rho_preleach_feed + rho_wash_water
        vol_flowrate_inflow = mass_flowrate_inflow/rho_std

        mass_flowrate_initial = std_underflow_mass_flowrate_initial + std_overflow_mass_flowrate_initial

        # Place holder for mass balance
        vol_flowrate_initial = mass_flowrate_initial/rho_std

        std_liq_volume_initial = self.single_tank_decantation_state_phase.get_value('liquid-volume', time)

        std_liq_volume = std_liq_volume_initial + \
                         (vol_flowrate_inflow - vol_flowrate_initial) * self.time_step

        if std_liq_volume > self.std_tank_volume:
            flow_residence_time = self.std_tank_volume/vol_flowrate_initial

            mass_flowrate = mass_flowrate_inflow + \
                            math.exp(-self.time_step/flow_residence_time) * \
                            (mass_flowrate_initial - mass_flowrate_inflow)
        else:
            mass_flowrate = 0.0

        u_over_o = self.std_underflow_overflow_mass_flowrate_ratio

        std_overflow_mass_flowrate = 1/(u_over_o + 1) * mass_flowrate
        std_underflow_mass_flowrate = u_over_o/(u_over_o + 1) * mass_flowrate

        # STD Math
        m_dot_pl_std = self.single_tank_decantation_preleach_feed_mass_flowrate
        c_pl_std = self.single_tank_decantation_preleach_feed_solids_massfrac
        m_dot_o = self.single_tank_decantation_raffinate_feed_mass_flowrate
        c_o = self.single_tank_decantation_raffinate_feed_solids_massfrac
        m_dot_u = std_underflow_mass_flowrate
        #c_u = std_underflow_solids_massfrac

        c_o = 100 + (c_pl_std-100)*2.78**(-1*time/1200)
        c_u = 9900 + (c_pl_std-9900)*2.78**(-1*time/1200)
        m_dot_o = (c_pl_std*m_dot_pl_std - c_u*m_dot_u)/c_o
        m_dot_u = (c_pl_std*m_dot_pl_std - c_o*m_dot_o)/c_u


        tmp_std_state = self.single_tank_decantation_state_phase.get_row(time)
        tmp_std_overflow = self.single_tank_decantation_overflow_phase.get_row(time)
        tmp_std_underflow = self.single_tank_decantation_underflow_phase.get_row(time)

        #---------------------
        # Evolve the CCD State
        #---------------------

        ccd_underflow_mass_flowrate_initial = self.ccd_underflow_phase.get_value('mass-flowrate',
                                                                                                     time)
        ccd_overflow_mass_flowrate_initial = self.ccd_overflow_phase.get_value('mass-flowrate', time)

        acidleach_feed_mass_flowrate = self.ccd_acidleach_feed_mass_flowrate
        rho_acidleach_feed = self.ccd_acidleach_feed_mass_density
        acidleach_feed_vol_flowrate = acidleach_feed_mass_flowrate / rho_acidleach_feed

        ccd_wash_water_vol_flowrate = self.wash_water_acidleach_feed_ratio * acidleach_feed_vol_flowrate
        rho_wash_water = self.wash_water_mass_density
        ccd_wash_water_mass_flowrate = ccd_wash_water_vol_flowrate * rho_wash_water
        wash_water_massfrac = 500 * unit.ppm

        # Ideal solution
        mass_flowrate_inflow = acidleach_feed_mass_flowrate + ccd_wash_water_mass_flowrate

        mass_flowrate_initial = ccd_underflow_mass_flowrate_initial + ccd_overflow_mass_flowrate_initial

        rho_std = rho_acidleach_feed + rho_wash_water

        # Place holder for mass balance
        vol_flowrate_initial = mass_flowrate_initial / rho_std

        flow_residence_time = self.ccd_tank_volume / vol_flowrate_initial

        mass_flowrate = mass_flowrate_inflow + \
                        math.exp(-self.time_step/flow_residence_time) * \
                       (mass_flowrate_initial - mass_flowrate_inflow)

        u_over_o = self.ccd_underflow_overflow_mass_flowrate_ratio

        ccd_overflow_mass_flowrate = 1 / (u_over_o + 1) * mass_flowrate
        ccd_underflow_mass_flowrate = u_over_o / (u_over_o + 1) * mass_flowrate

        # CCD Math
        m_dot_al = self.ccd_acidleach_feed_mass_flowrate
        c_al = self.ccd_acidleach_feed_solids_massfrac
        m_dot_t = self.ccd_underflow_mass_flowrate
        c_t = self.ccd_underflow_solids_massfrac
        m_dot_pl_ccd = self.ccd_overflow_mass_flowrate
        c_pl_ccd = self.ccd_overflow_solids_massfrac

        c_t = 9801 + (c_al - 9801)*2.78**(-1*time/1200)
        c_pl_ccd = 99 + (c_al - 99)*2.78**(-1*time/1200)
        m_dot_pl_ccd = (c_al*m_dot_al + 500*(1.69/3.28)*m_dot_al - c_t*m_dot_t)/c_pl_ccd
        m_dot_t = (c_al*m_dot_al + 500*(1.69/3.28)*m_dot_al - c_pl_ccd*m_dot_pl_ccd)/c_t

        tmp_ccd_overflow = self.ccd_overflow_phase.get_row(time)
        tmp_ccd_underflow = self.ccd_underflow_phase.get_row(time)

        #----------------------------
        # Evolve the Filtration State
        #----------------------------

        filtration_filtrate_mass_flowrate_initial = self.filtration_filtrate_phase.get_value('mass-flowrate',
                                                                                 time)
        filtration_slurry_mass_flowrate_initial = self.filtration_slurry_phase.get_value('mass-flowrate', time)

        filtration_feed_mass_flowrate = std_overflow_mass_flowrate
        rho_filtration_feed = rho_std
        filtration_feed_vol_flowrate = filtration_feed_mass_flowrate / rho_filtration_feed

        # Ideal solution
        mass_flowrate_inflow = filtration_feed_mass_flowrate

        mass_flowrate_initial = filtration_filtrate_mass_flowrate_initial + filtration_slurry_mass_flowrate_initial

        # Filtration Math
        m_dot_sl = self.filtration_slurry_mass_flowrate
        c_sl = self.filtration_slurry_solids_massfrac
        m_dot_f = self.filtration_filtrate_mass_flowrate
        c_f = self.filtration_filtrate_solids_massfrac

        c_f = 10 * unit.ppm
        c_sl = 990 * unit.ppm
        m_dot_f = 0.991*m_dot_o
        m_dot_sl = 0.009*m_dot_o

        tmp_filtrate = self.filtration_filtrate_phase.get_row(time)
        tmp_slurry = self.filtration_slurry_phase.get_row(time)

        time += self.time_step

        self.single_tank_decantation_state_phase.add_row(time, tmp_std_state)
        self.single_tank_decantation_state_phase.set_value('liquid-volume', std_liq_volume, time)

        self.single_tank_decantation_overflow_phase.add_row(time, tmp_std_overflow)
        self.single_tank_decantation_overflow_phase.set_value('mass-flowrate',
                                                              std_overflow_mass_flowrate,
                                                              time)

        self.single_tank_decantation_underflow_phase.add_row(time, tmp_std_underflow)
        self.single_tank_decantation_underflow_phase.set_value('mass-flowrate',
                                                               std_underflow_mass_flowrate, time)

        self.ccd_overflow_phase.add_row(time, tmp_ccd_overflow)
        self.ccd_overflow_phase.set_value('mass-flowrate', ccd_overflow_mass_flowrate, time)

        self.ccd_underflow_phase.add_row(time, tmp_ccd_underflow)
        self.ccd_underflow_phase.set_value('mass-flowrate', ccd_underflow_mass_flowrate, time)

        self.filtration_filtrate_phase.add_row(time, tmp_filtrate)
        self.filtration_filtrate_phase.set_value('mass-flowrate', m_dot_f, time)

        self.filtration_slurry_phase.add_row(time, tmp_slurry)
        self.filtration_slurry_phase.set_value('mass-flowrate', m_dot_sl, time)


        return time
