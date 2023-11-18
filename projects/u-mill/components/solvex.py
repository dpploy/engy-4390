#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""
Cortix Module
This module is a model of the Solvent Extraction process in the White Mesa Uranium Milling Plant
"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

import unit

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class Solvex(Module):
    """Steam generator.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: Filtration, Precipitation.
    See instance attribute `port_names_expected`.

    """
    
    def __init__(self):
        #Original: def __init__(self, primary_inflow_temp=20+273.15, secondary_inflow_temp=20+273.15):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['primary-inflow', 'primary-outflow',
                                    'secondary-inflow', 'secondary-outflow']

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
        '''
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
        '''
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
        
        # Derived quantities
        '''
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

        # Primary outflow phase history
        quantities = list()

        '''
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
            self.__call_ports(time)

        self.end_time = time # correct the final time if needed

    def __call_ports(self, time):

        # Interactions in the uranium-inflow port
        #----------------------------------------
        # One way "from" uranium-inflow

        # Receive from
        if self.get_port('uranium-inflow').connected_port:

            self.send(time, 'uranium-inflow')

            (check_time, primary_inflow) = self.recv('uranium-inflow')
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
