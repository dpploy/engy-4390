#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module.
   Decantation/Filtration process in the White Mesa Uranium Milling Plant.
   Add info here:...

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

        self.feed_aqueous_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.feed_solid_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.raffinate_aqueous_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.raffinate_solid_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.tailings_solid_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.tailings_aqueous_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_raffinate_aqueous_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_raffinate_solid_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_tailings_aqueous_mass_flowrate = 1.0 * unit.liter/unit.minute
        self.filtration_tailings_solid_mass_flowrate = 1.0 * unit.liter/unit.minute
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

        # Aqueous feed decantation history
        quantities = list()
        species = list()

        feed_aqueous_flowrate = Quantity(name='flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.feed_aqueous_mass_flowrate,
                        latex_name=r'$\dot{m}_1$',
                        info='Clarification Feed Aqueous Mass Flowrate')
        quantities.append(feed_aqueous_flowrate)

        uo2so434minus_feed = Species( name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(feed)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)
        
        h2o_feed = Species(name='H2O',formula_name='H2O(feed)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        h2so4_feed = Species(name='H2SO4',formula_name='H2SO4(feed)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_feed)
        self.feed_aqueous_phase = Phase(time_stamp=self.initial_time,
                                        time_unit='s', quantities=quantities, species=species)
        # Solid feed decantation history (Iron, Gold, and Copper impurities)
        quantities = list()
        species = list()#Is it proper to rest the species list too?

        feed_solid_flowrate = Quantity(name='flowrate',
                        formal_name='mdot', unit='kg/s',
                        value=self.feed_solid_mass_flowrate,
                        latex_name=r'$\dot{m}_2$',
                        info='Clarification Feed Solid Mass Flowrate')
        quantities.append(feed_solid_flowrate)

        iron_feed = Species(name='Fe', formula_name='Fe(feed)',
                           atoms=['Fe'],
                           info='Fe')
        species.append(iron_feed)

        copper_feed = Species(name='Cu', formula_name='Cu(feed)',
                           atoms=['Cu'],
                           info='Cu')
        species.append(copper_feed)

        gold_feed = Species(name='Au', formula_name='Au(feed)',
                           atoms=['Au'],
                           info='Au')
        species.append(gold_feed)
        self.feed_solid_phase = Phase(time_stamp=self.initial_time,
                                       time_unit='s', quantities=quantities, species=species)

        # Aqueous raffinate decantation history
        quantities = list()
        species = list()

        raffinate_aqueous_flowrate = Quantity(name='flowrate',
                                         formal_name='mdot', unit='kg/s',
                                         value=self.raffinate_aqueous_mass_flowrate,
                                         latex_name=r'$\dot{m}_3$',
                                         info='Clarification Raffinate Aqueous Mass Flowrate')
        quantities.append(raffinate_aqueous_flowrate)

        uo2so434minus_raffinate = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(raffinate)',
                                     atoms=['U', '2*O', '3*S', '12*O'],
                                     info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_raffinate)

        h2o_raffinate = Species(name='H2O', formula_name='H2O(raffinate)',
                           atoms=['2*H', 'O'],
                           info='H2O')
        species.append(h2o_raffinate)

        h2so4_raffinate = Species(name='H2SO4', formula_name='H2SO4(raffinate)',
                             atoms=['2*H', 'S', '4*O'],
                             info='H2SO4')
        species.append(h2so4_raffinate)
        self.raffinate_aqueous_phase = Phase(time_stamp=self.initial_time,
                                        time_unit='s', quantities=quantities, species=species)

        # Solid raffinate decantation history
        quantities = list()
        species = list()  # Is it proper to rest the species list too?

        raffinate_solid_flowrate = Quantity(name='flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=self.raffinate_solid_mass_flowrate,
                                       latex_name=r'$\dot{m}_4$',
                                       info='Clarification raffinate Solid Mass Flowrate')
        quantities.append(raffinate_solid_flowrate)

        iron_raffinate = Species(name='Fe', formula_name='Fe(raffinate)',
                            atoms=['Fe'],
                            info='Fe')
        species.append(iron_raffinate)

        copper_raffinate = Species(name='Cu', formula_name='Cu(raffinate)',
                              atoms=['Cu'],
                              info='Cu')
        species.append(copper_raffinate)

        gold_raffinate = Species(name='Au', formula_name='Au(raffinate)',
                            atoms=['Au'],
                            info='Au')
        species.append(gold_raffinate)
        self.raffinate_solid_phase = Phase(time_stamp=self.initial_time,
                                      time_unit='s', quantities=quantities, species=species)

        # Solid tailings decantation history
        quantities = list()
        species = list()  # Is it proper to rest the species list too?

        tailings_solid_flowrate = Quantity(name='flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=self.tailings_solid_mass_flowrate,
                                       latex_name=r'$\dot{m}_4$',
                                       info='Clarification tailings Solid Mass Flowrate')
        quantities.append(tailings_solid_flowrate)

        iron_tailings = Species(name='Fe', formula_name='Fe(tailings)',
                            atoms=['Fe'],
                            info='Fe')
        species.append(iron_tailings)

        copper_tailings = Species(name='Cu', formula_name='Cu(tailings)',
                              atoms=['Cu'],
                              info='Cu')
        species.append(copper_tailings)

        gold_tailings = Species(name='Au', formula_name='Au(tailings)',
                            atoms=['Au'],
                            info='Au')
        species.append(gold_tailings)
        self.tailings_solid_phase = Phase(time_stamp=self.initial_time,
                                      time_unit='s', quantities=quantities, species=species)

        # Aqueous tailings decantation history
        quantities = list()
        species = list()

        tailings_aqueous_flowrate = Quantity(name='flowrate',
                                              formal_name='mdot', unit='kg/s',
                                              value=self.tailings_aqueous_mass_flowrate,
                                              latex_name=r'$\dot{m}_3$',
                                              info='Clarification tailings Aqueous Mass Flowrate')
        quantities.append(raffinate_aqueous_flowrate)

        uo2so434minus_tailings = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(raffinate)',
                                          atoms=['U', '2*O', '3*S', '12*O'],
                                          info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_tailings)

        h2o_tailings = Species(name='H2O', formula_name='H2O(tailings)',
                                atoms=['2*H', 'O'],
                                info='H2O')
        species.append(h2o_tailings)

        h2so4_tailings = Species(name='H2SO4', formula_name='H2SO4(tailings)',
                                  atoms=['2*H', 'S', '4*O'],
                                  info='H2SO4')
        species.append(h2so4_tailings)
        self.tailings_aqueous_phase = Phase(time_stamp=self.initial_time,
                                             time_unit='s', quantities=quantities, species=species)

        # Aqueous feed filtration history
        # Same as Decantation aqueous raffinate

        # Solid feed filtration history
        # Same as Decantation solid raffinate

        # Aqueous raffinate filtration history
        quantities = list()
        species = list()

        filtration_raffinate_aqueous_flowrate = Quantity(name='flowrate',
                                              formal_name='mdot', unit='kg/s',
                                              value=self.filtration_raffinate_aqueous_mass_flowrate,
                                              latex_name=r'$\dot{m}_3$',
                                              info='Filtration Raffinate Aqueous Mass Flowrate')
        quantities.append(filtration_raffinate_aqueous_flowrate)

        uo2so434minus_filtration_raffinate = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(filtration raffinate)',
                                          atoms=['U', '2*O', '3*S', '12*O'],
                                          info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_filtration_raffinate)

        h2o_filtration_raffinate = Species(name='H2O', formula_name='H2O(filtration raffinate)',
                                atoms=['2*H', 'O'],
                                info='H2O')
        species.append(h2o_filtration_raffinate)

        h2so4_filtration_raffinate = Species(name='H2SO4', formula_name='H2SO4(filtration raffinate)',
                                  atoms=['2*H', 'S', '4*O'],
                                  info='H2SO4')
        species.append(h2so4_filtration_raffinate)
        self.filtration_raffinate_aqueous_phase = Phase(time_stamp=self.initial_time,
                                             time_unit='s', quantities=quantities, species=species)

        # Solid raffinate filtration history
        # Essentially 0% wt

        # Aqueous tailings filtration history
        quantities = list()
        species = list()

        filtration_tailings_aqueous_flowrate = Quantity(name='flowrate',
                                                         formal_name='mdot', unit='kg/s',
                                                         value=self.filtration_tailings_aqueous_mass_flowrate,
                                                         latex_name=r'$\dot{m}_3$',
                                                         info='Filtration tailings Aqueous Mass Flowrate')
        quantities.append(filtration_tailings_aqueous_flowrate)

        uo2so434minus_filtration_tailings = Species(name='UO2-(SO4)3^4-',
                                                     formula_name='UO2(SO4)3^4-(filtration tailings)',
                                                     atoms=['U', '2*O', '3*S', '12*O'],
                                                     info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_filtration_tailings)

        h2o_filtration_tailings = Species(name='H2O', formula_name='H2O(filtration tailings)',
                                           atoms=['2*H', 'O'],
                                           info='H2O')
        species.append(h2o_filtration_tailings)

        h2so4_filtration_tailings = Species(name='H2SO4', formula_name='H2SO4(filtration tailings)',
                                             atoms=['2*H', 'S', '4*O'],
                                             info='H2SO4')
        species.append(h2so4_filtration_tailings)
        self.filtration_tailings_aqueous_phase = Phase(time_stamp=self.initial_time,
                                                        time_unit='s', quantities=quantities, species=species)

        # Solid tailings filtration history
        quantities = list()
        species = list()  # Is it proper to rest the species list too?

        filtration_tailings_solid_flowrate = Quantity(name='flowrate',
                                           formal_name='mdot', unit='kg/s',
                                           value=self.filtration_tailings_solid_mass_flowrate,
                                           latex_name=r'$\dot{m}_4$',
                                           info='Filtration tailings Solid Mass Flowrate')
        quantities.append(filtration_tailings_solid_flowrate)

        iron_filtration_tailings = Species(name='Fe', formula_name='Fe(filtration tailings)',
                                atoms=['Fe'],
                                info='Fe')
        species.append(iron_filtration_tailings)

        copper_filtration_tailings = Species(name='Cu', formula_name='Cu(filtration tailings)',
                                  atoms=['Cu'],
                                  info='Cu')
        species.append(copper_filtration_tailings)

        gold_filtration_tailings = Species(name='Au', formula_name='Au(filtration tailings)',
                                atoms=['Au'],
                                info='Au')
        species.append(gold_filtration_tailings)
        self.filtration_tailings_solid_phase = Phase(time_stamp=self.initial_time,
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
