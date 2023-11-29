#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module.
   Decantation/Filtration process in the White Mesa Uranium Milling Plant.
   Add info here:...

   + Decantation:
       - feed concentration of solids:       100000 ppm
       - wash water concentration of solids: 500 ppm
       - overflow concentration of solids:   100 ppm
       - # of thickeners:                    4
       - underflow concentration of solids:  99000 ppm


   + Wash water (volume? mass?) flowrate:
   + Feed (aqueous/solids) (volume? mass?) flow rate:

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

        self.port_names_expected = ['decantation-inflow', 'decantation-outflow',
                                     'drum-inflow', 'drum-outflow'] #Need to define these

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
        self.wash_water_flowrate = 1.0 * unit.liter/unit.minute

        '''
        Example

        self.discard_tau_recording_before = 2*unit.minute
        self.heat_transfer_area = 1665.57*unit.meter**2

        self.helicoil_outer_radius = 16/2*unit.milli*unit.meter
        self.helicoil_tube_wall = 0.9*unit.milli*unit.meter
        self.helicoil_inner_radius = self.helicoil_outer_radius - self.helicoil_tube_wall
        self.helicoil_length = 22.3*unit.meter
        self.n_helicoil_tubes = 1380

        self.wall_temp_delta_primary = 1.5*unit.K
        self.wall_temp_delta_secondary = 1.5*unit.K

        self.iconel690_k = 12.1*unit.watt/unit.meter/unit.kelvin

        self.helix_to_cylinder = 1./.928

        self.secondary_volume = math.pi * self.helicoil_inner_radius**2 * \
                                self.helicoil_length * self.n_helicoil_tubes *\
                                self.helix_to_cylinder

        self.primary_volume = 0.5 * self.secondary_volume

        # Ratio of the tube bundle pithc transverse to flow to parallel to flow
        self.tube_bundle_pitch_ratio = 1.5  # st/sl
        '''

        # Initialization

        # Decantation          
        self.decantation_feed_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.decantation_feed_solids_conc = 100 * unit.ppm

        self.decantation_underflow_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.decantation_underflow_solids_conc = 100 * unit.ppm

        self.decantation_overflow_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.decantation_overflow_solids_conc = 100 * unit.ppm

        # Filtration

        self.filtration_underflow_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_underflow_solids_conc = 100 * unit.ppm



        '''
        Example

        self.primary_inflow_temp = primary_inflow_temp

        self.primary_pressure = 127.6*unit.bar

        self.primary_mass_flowrate = 4.66e6*unit.lb/unit.hour

        self.primary_outflow_temp = self.primary_inflow_temp #- 2*unit.K

        self.secondary_inflow_temp = secondary_inflow_temp

        self.secondary_pressure = 34*unit.bar
        self.secondary_mass_flowrate = 67*unit.kg/unit.second

        self.secondary_outflow_temp = self.secondary_inflow_temp #- 2*unit.K

        self.secondary_outflow_quality = 0 # running value of quality
        '''
        '''
        # Derived quantities
        self.rho_p = 0.0
        self.cp_p = 0.0
        self.mu_p = 0.0
        self.k_p = 0.0
        self.prtl_p = 0.0
        self.rey_p = 0.0
        self.nusselt_p = 0.0

        self.rho_s = 0.0
        self.cp_s = 0.0
        self.mu_s = 0.0
        self.k_s = 0.0
        self.prtl_s = 0.0
        self.rey_s = 0.0
        self.nusselt_s = 0.0

        self.heat_sink_pwr = 0.0
        '''

        #***************************************************************************************
        # D E C A N T A T I O N 
        #***************************************************************************************

        # Decantation Feed Phase History (leaching outflow)
        quantities = list()
        species = list()

        feed_mass_flowrate = Quantity(name='mass_flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.decantation_feed_mass_flowrate,
                        latex_name=r'$\dot{m}_1$',
                        info='Decantation Feed Mass Flowrate')
        quantities.append(feed_mass_flowrate)

        feed_solids_conc = Quantity(name='solids_conc',
                        formal_name='solids_conc', unit='ppm',
                        value=self.decantation_feed_solids_conc,
                        latex_name=r'$C_1$',
                        info='Decantation Feed Solids Concentration')

        quantities.append(feed_solids_conc)

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

        self.feed_phase = Phase(time_stamp=self.initial_time,
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

        underflow_solids_conc = Quantity(name='solids_conc',
                        formal_name='solids_conc', unit='ppm',
                        value=self.decantation_underflow_solids_conc,
                        latex_name=r'$C_1$',
                        info='Decantation Underflow Solids Concentration')

        quantities.append(underflow_solids_conc)

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

        self.underflow_phase = Phase(time_stamp=self.initial_time,
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

        overflow_solids_conc = Quantity(name='solids_conc',
                        formal_name='solids_conc', unit='ppm',
                        value=self.decantation_overflow_solids_conc,
                        latex_name=r'$C_1$',
                        info='Decantation Overflow Solids Concentration')

        quantities.append(overflow_solids_conc)

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

        self.overflow_phase = Phase(time_stamp=self.initial_time,
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

        slurry_solids_conc = Quantity(name='solids_conc',
                        formal_name='solids_conc', unit='ppm',
                        value=self.filtration_slurry_solids_conc,
                        latex_name=r'$C_1$',
                        info='Filtration Slurry Solids Concentration')

        quantities.append(slurry_solids_conc)

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

        self.slurry_phase = Phase(time_stamp=self.initial_time,
                                  time_unit='s', quantities=quantities, species=species)





        '''
        quantities = list()
        temp = Quantity(name='temp',
                        formal_name='T_1', unit='K',
                        value=self.primary_outflow_temp,
                        latex_name=r'$T_1$',
                        info='Steamer Primary Outflow Temperature')

        quantities.append(temp)

        flowrate = Quantity(name='flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.primary_mass_flowrate,
                        latex_name=r'$\dot{m}_1$',
                        info='Steamer Primary Mass Flowrate')

        quantities.append(flowrate)

        self.primary_outflow_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s', quantities=quantities)

        # Secondary inflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='m2i', unit='kg/s',
                            value=self.secondary_mass_flowrate,
                            latex_name=r'$\dot{m}_{2,in}$',
                            info='Steamer Secondary Inflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T2i', unit='K',
                        value=self.secondary_inflow_temp,
                        latex_name=r'$T_{2,in}$',
                        info='Steamer Secondary Inflow Temperature')

        quantities.append(temp)

        self.secondary_inflow_phase = Phase(time_stamp=self.initial_time,
                                            time_unit='s', quantities=quantities)

        # Secondary outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='m_2', unit='kg/s',
                            value=self.secondary_mass_flowrate,
                            latex_name=r'$\dot{m}_2$',
                            info='Steamer Secondary Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.secondary_outflow_temp,
                        latex_name=r'$T_2$',
                        info='Steamer Secondary Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='Pa',
                         value=self.secondary_pressure,
                         latex_name=r'$P_2$',
                         info='Steamer Secondary Outflow Pressure')

        quantities.append(press)

        quality = Quantity(name='quality',
                         formal_name='X', unit='',
                         value=self.secondary_outflow_quality,
                         latex_name=r'$\chi$',
                         info='Steamer Secondary Outflow Quality')

        quantities.append(quality)

        self.secondary_outflow_phase = Phase(time_stamp=self.initial_time,
                                             time_unit='s', quantities=quantities)

        # State phase history
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
            #self.__call_ports(time)

        self.end_time = time # correct the final time if needed



    '''
    def __call_ports(self, time):

        # Interactions in the decantation-inflow port
        #----------------------------------------
        # One way "from" decantation-inflow

        # Receive from
        if self.get_port('decantation-inflow').connected_port:

            self.send(time, 'decantation-inflow')

            (check_time, primary_inflow) = self.recv('decantation-inflow')
            assert abs(check_time-time) <= 1e-6

            self.primary_inflow_temp = primary_inflow['temperature']
            self.primary_ressure = primary_inflow['pressure']
            self.primary_mass_flowrate = primary_inflow['mass_flowrate']

        # Interactions in the secondary-inflow port
        #----------------------------------------
        # One way "from" secondary-inflow

        # Receive from
        if self.get_port('secondary-inflow').connected_port:

            self.send(time, 'secondary-inflow')

            (check_time, secondary_inflow) = self.recv('secondary-inflow')
            assert abs(check_time-time) <= 1e-6

            self.secondary_inflow_temp = secondary_inflow['temperature']
            self.secondary_pressure = secondary_inflow['pressure']
            self.secondary_mass_flowrate = secondary_inflow['mass_flowrate']

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # One way "to" primary-outflow

        # Send to
        if self.get_port('primary-outflow').connected_port:

            msg_time = self.recv('primary-outflow')

            temp = self.primary_outflow_phase.get_value('temp', msg_time)

            primary_outflow = dict()
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_pressure
            primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
            primary_outflow['quality'] = 0.0

            self.send((msg_time, primary_outflow), 'primary-outflow')

        # Interactions in the secondary-outflow port
        #-----------------------------------------
        # One way "to" secondary-outflow

        # Send to
        if self.get_port('secondary-outflow').connected_port:

            msg_time = self.recv('secondary-outflow')

            temp = self.secondary_outflow_phase.get_value('temp', msg_time)
            press = self.secondary_outflow_phase.get_value('pressure', msg_time)
            flowrate = self.secondary_outflow_phase.get_value('flowrate', msg_time)

            secondary_outflow = dict()
            secondary_outflow['temperature'] = temp
            secondary_outflow['pressure'] = press
            secondary_outflow['mass_flowrate'] = flowrate
            secondary_outflow['total_heat_power'] = -self.heat_sink_pwr

            self.send((msg_time, secondary_outflow), 'secondary-outflow')
    '''

    def __call_ports(self, time):

        # Interactions in the decantation-inflow port
        #----------------------------------------
        # One way "from" decantation-inflow

        # Receive from
        if self.get_port('decantation-inflow').connected_port:

            self.send(time, 'decantation-inflow')

            (check_time, primary_inflow) = self.recv('decantation-inflow')
            assert abs(check_time-time) <= 1e-6

            self.primary_inflow_temp = primary_inflow['temperature']
            self.primary_ressure = primary_inflow['pressure']
            self.primary_mass_flowrate = primary_inflow['mass_flowrate']

        # Interactions in the secondary-inflow port
        #----------------------------------------
        # One way "from" secondary-inflow

        # Receive from
        if self.get_port('secondary-inflow').connected_port:

            self.send(time, 'secondary-inflow')

            (check_time, secondary_inflow) = self.recv('secondary-inflow')
            assert abs(check_time-time) <= 1e-6

            self.secondary_inflow_temp = secondary_inflow['temperature']
            self.secondary_pressure = secondary_inflow['pressure']
            self.secondary_mass_flowrate = secondary_inflow['mass_flowrate']

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # One way "to" primary-outflow

        # Send to
        if self.get_port('primary-outflow').connected_port:

            msg_time = self.recv('primary-outflow')

            temp = self.primary_outflow_phase.get_value('temp', msg_time)

            primary_outflow = dict()
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_pressure
            primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
            primary_outflow['quality'] = 0.0

            self.send((msg_time, primary_outflow), 'primary-outflow')

        # Interactions in the secondary-outflow port
        #-----------------------------------------
        # One way "to" secondary-outflow

        # Send to
        if self.get_port('secondary-outflow').connected_port:

            msg_time = self.recv('secondary-outflow')

            temp = self.secondary_outflow_phase.get_value('temp', msg_time)
            press = self.secondary_outflow_phase.get_value('pressure', msg_time)
            flowrate = self.secondary_outflow_phase.get_value('flowrate', msg_time)

            secondary_outflow = dict()
            secondary_outflow['temperature'] = temp
            secondary_outflow['pressure'] = press
            secondary_outflow['mass_flowrate'] = flowrate
            secondary_outflow['total_heat_power'] = -self.heat_sink_pwr

            self.send((msg_time, secondary_outflow), 'secondary-outflow')
