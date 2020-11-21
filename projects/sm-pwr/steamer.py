#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

import unit

import iapws.iapws97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class Steamer(Module):
    """Steam generator.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: reactor, turbine.
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
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

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters
        #self.ntu = 3.76
        self.effectiveness = 0.949
        self.c_min = 1110.5 #kJ/K
        self.c_hot = 3729.7 #kJ/K
        self.real_cp_secondary =   4.2980 #kJ/kg-K
        self.t_sat = 516 #K
        self.q_vap = 145350.6 #kj
        self.mw_water = 18.02 #g/mol

        # Initialization
        self.primary_inflow_pressure = 12.8*unit.mega*unit.pascal
        self.primary_inflow_temp = 30+273.15
        self.primary_inflow_mass_flowrate = 666

        self.secondary_inflow_pressure = 3.4*unit.mega*unit.pascal
        self.secondary_inflow_temp = 20+273.15
        self.secondary_inflow_mass_flowrate = 67

        self.primary_outflow_temp = 20 + 273.15
        self.primary_outflow_mass_flowrate = 67
        self.primary_outflow_pressure = 12.8*unit.mega*unit.pascal

        self.secondary_outflow_mass_flowrate = 67
        self.secondary_outflow_temp = 20 + 273.15
        self.secondary_outflow_pressure = 3.4*unit.mega*unit.pascal

        # Primary outflow phase history
        quantities = list()

        temp = Quantity(name='temp',
                        formal_name='T_1', unit='K',
                        value=self.primary_outflow_temp,
                        latex_name=r'$T_1$',
                        info='Steamer Primary Outflow Temperature')

        quantities.append(temp)

        self.primary_outflow_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s', quantities=quantities)

        # Secondary outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.secondary_outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Steamer Secondary Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.secondary_outflow_temp,
                        latex_name=r'$T_2$',
                        info='Steamer Secondary Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='p_out', unit='Pa',
                         value=self.secondary_outflow_pressure,
                         latex_name=r'$p_out$',
                         info='Steamer Secondary Outflow Pressure')

        quantities.append(press)

        self.secondary_outflow_phase = Phase(time_stamp=self.initial_time,
                                             time_unit='s', quantities=quantities)

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

            # Communicate information
            #------------------------
            self.__call_ports(time)

            # Evolve one time step
            #---------------------

            time = self.__step(time)

    def __call_ports(self, time):

        # Interactions in the primary-inflow port
        #----------------------------------------
        # One way "from" primary-inflow

        # Receive from
        if self.get_port('primary-inflow').connected_port:

            self.send(time, 'primary-inflow')

            (check_time, primary_inflow) = self.recv('primary-inflow')
            assert abs(check_time-time) <= 1e-6

            primary_inflow['temperature'] = self.primary_inflow_temp
            primary_inflow['pressure'] = self.primary_inflow_pressure
            primary_inflow['mass_flowrate'] = self.primary_inflow_mass_flowrate

        # Interactions in the secondary-inflow port
        #----------------------------------------
        # One way "from" secondary-inflow

        # Receive from
        if self.get_port('secondary-inflow').connected_port:

            self.send(time, 'secondary-inflow')

            (check_time, secondary_inflow) = self.recv('secondary-inflow')
            assert abs(check_time-time) <= 1e-6

            self.secondary_inflow_temp = secondary_inflow['temperature']
            self.secondary_inflow_pressure = secondary_inflow['pressure']
            self.secondary_inflow_mass_flowrate = secondary_inflow['temperature']

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # One way "to" primary-outflow

        # Send to
        if self.get_port('primary-outflow').connected_port:

            msg_time = self.recv('primary-outflow')

            temp = self.primary_outflow_phase.get_value('temp', msg_time)
            primary_outflow = dict()
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_outflow_pressure
            primary_outflow['flowrate'] = self.primary_outflow_mass_flowrate
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
            secondary_outflow['flowrate'] = flowrate

            self.send((msg_time, secondary_outflow), 'secondary-outflow')

    def __step(self, time=0.0):

        p_in = self.primary_inflow_pressure
        p_out = self.secondary_outflow_pressure
        t_h_i = self.primary_inflow_temp
        t_c_i = self.secondary_inflow_temp
        c_min = self.c_min
        cp_out = self.real_cp_secondary
        # Calculations

        q = self.effectiveness*c_min*(t_h_i-t_c_i)

        t_h_o = t_h_i - (q/self.c_hot)
        t_c_o = t_c_i + (q/c_min)

        q_s = c_min*(t_c_o - t_c_i)
        q_req = cp_out*(self.t_sat-t_c_i)

        # If no vaporization
        if q_s < q_req: #kW

            t_exit = q_s/(cp_out*self.secondary_inflow_mass_flowrate) + t_c_i

        # If some vaporization
        elif q_s < self.q_vap:

            t_exit = self.t_sat

        # If complete vaporization - calculate gas tenperautre with cp function
        else:

            q_heat_steam = q_s - self.q_vap
            delta_h = q_heat_steam/self.secondary_inflow_mass_flowrate

            a = 3.470
            b =  1.45*10e-3
            c = (delta_h*self.mw_water/unit.R) + a*self.t_sat + b*(self.t_sat**2)
            coeff = [b,a,-c]
            roots = np.roots(coeff)

            if roots[0] > 0:
                t_exit = roots[0]
            else:
                t_exit = roots[1]

        # Update phases
        primary_outflow = self.primary_outflow_phase.get_row(time)
        secondary_outflow = self.secondary_outflow_phase.get_row(time)

        time += self.time_step

        self.primary_outflow_phase.add_row(time, primary_outflow)
        self.primary_outflow_phase.set_value('temp', t_h_o, time)

        self.secondary_outflow_phase.add_row(time, secondary_outflow)
        self.secondary_outflow_phase.set_value('temp', t_exit, time)
        self.secondary_outflow_phase.set_value('flowrate', self.secondary_inflow_mass_flowrate, time)
        self.secondary_outflow_phase.set_value('pressure', p_out, time)

        return time
