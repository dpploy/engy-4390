#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import unit

import iapws.iapws97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class Turbine(Module):
    """Turbine.

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

        self.port_names_expected = ['inflow', 'outflow']

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters

        # Initialization
        self.turbine_efficiency = 0.7784

        self.inflow_pressure = 1.0*unit.bar
        self.inflow_temp = 20+273.15
        self.inflow_mass_flowrate = 0.0*unit.kg/unit.second

        self.outflow_temp = 20+272.15
        self.outflow_mass_flowrate = 0.0
        self.outflow_pressure = 0.008066866*unit.mega*unit.pascal
        self.outflow_quality = 0.0

        # Outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q$',
                            info='Turbine Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T$',
                        info='Turbine Outflow Temperature')

        quantities.append(temp)

        temp = Quantity(name='pressure',
                        formal_name='P', unit='Pa',
                        value=self.outflow_pressure,
                        latex_name=r'$P$',
                        info='Turbine Outflow Pressure')

        quantities.append(temp)

        quality = Quantity(name='quality',
                         formal_name='X', unit='',
                         value=self.outflow_quality,
                         latex_name=r'$X$',
                         info='Turbine Exit Steam Quality')

        quantities.append(quality)

        self.outflow_phase = Phase(time_stamp=self.initial_time,
                                   time_unit='s', quantities=quantities)

        # Turbine phase history
        quantities = list()

        power = Quantity(name='power',
                         formal_name='W_s', unit='W_e',
                         value=0.0,
                         latex_name=r'$W_s$',
                         info='Turbine Power')

        quantities.append(power)

        self.state_phase = Phase(time_stamp=self.initial_time,
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

        # Interactions in the outflow port
        #-----------------------------------------
        # one way "to" outflow

        # send to
        if self.get_port('outflow').connected_port:

            msg_time = self.recv('outflow')

            temp = self.outflow_phase.get_value('temp', msg_time)
            pressure = self.outflow_phase.get_value('pressure', msg_time)
            outflow = dict()
            outflow['temperature'] = temp
            outflow['pressure'] = pressure
            outflow['flowrate'] = self.outflow_mass_flowrate

            self.send((msg_time, outflow), 'outflow')

        # Interactions in the inflow port
        #----------------------------------------
        # one way "from" inflow

        # receive from
        if self.get_port('inflow').connected_port:

            self.send(time, 'inflow')

            (check_time, inflow) = self.recv('inflow')
            assert abs(check_time-time) <= 1e-6

            self.inflow_temp = inflow['temperature']
            self.inflow_pressure = inflow['pressure']
            self.inflow_mass_flowrate = inflow['mass_flowrate']

    def __step(self, time=0.0):

        # Check for valid intervals

        p_in_min = 6.1121e-4*unit.mega*unit.pascal
        p_in_max = 22.064*unit.mega*unit.pascal
        assert p_in_min <= self.inflow_pressure <= p_in_max

        assert 20+273.15 <= self.inflow_temp <= 800+273.15

        # Get state values
        p_in_MPa = self.inflow_pressure/unit.mega/unit.pascal
        temp_in = self.inflow_temp
        p_out_MPa = self.outflow_pressure/unit.mega/unit.pascal
        m_dot_2 = self.inflow_mass_flowrate

        #############################################

        #if entering stream is not steam (valve closed scenario)
        if temp_in < steam_table._TSat_P(p_in_MPa):
            t_runoff = steam_table._TSat_P(p_in_MPa)
            power = 0
            quality = 0

        else:
            s_2_prime = steam_table._Region2(temp_in, p_in_MPa)['s']
            h_1 = steam_table._Region2(temp_in, p_in_MPa)['h']
            bubl = steam_table._Region4(p_out_MPa, 0)
            dew = steam_table._Region4(p_out_MPa, 1)
            bubl_entropy = bubl['s']
            dew_entropy = dew['s']
            bubl_enthalpy = bubl['h']
            dew_enthalpy = dew['h']

            #if the ideal runoff is two-phase mixture:
            if bubl_entropy < s_2_prime < dew_entropy:

                quality = (s_2_prime - bubl_entropy) / (dew_entropy - bubl_entropy)

                h_2_prime = bubl_enthalpy + quality * (dew_enthalpy - bubl_enthalpy)

            #if ideal run off is superheated
            elif s_2_prime > dew_entropy:

                t_ideal = steam_table._Backward2_T_Ps(p_out_MPa, s_2_prime)

                h_2_prime = steam_table._Region2(t_ideal, p_out_MPa)['h']

                quality = 1

            #else ideal run off is subcooled
            else:
                t_ideal = steam_table._Backward1_T_Ps(p_out_MPa, s_2_prime)

                h_2_prime = steam_table._Region1(t_ideal, p_out_MPa)['h']

                quality = 0


            #calculate the real runoff enthalpy
            w_ideal = h_1 - h_2_prime  #on a per mass basis
            assert w_ideal > 0
            w_real = w_ideal * self.turbine_efficiency
            h_real = h_1 - w_ideal
            assert h_real > 0
            if w_real < 0:
                w_real = 0

            #if run off is actually subcooled
            if h_real < bubl_enthalpy:

                t_runoff = steam_table._Backward1_T_Ph(p_out_MPa, h_real)

                quality = 0  # subcooled liquid

            #if run off is actually superheated
            elif h_real > dew_enthalpy:

                t_runoff = steam_table._Backward2_T_Ph(p_out_MPa, h_real)

                quality = 1 # superheated steam
            #else run off is actually in two phase region
            else:

                quality = (h_real - bubl_enthalpy)/(dew_enthalpy - bubl_enthalpy)

                t_runoff = steam_table._Region4(p_out_MPa, quality)['T']

            power = m_dot_2 * w_real

        #update state variables
        turbine_outflow = self.outflow_phase.get_row(time)
        turbine = self.state_phase.get_row(time)

        time += self.time_step

        self.outflow_phase.add_row(time, turbine_outflow)

        self.outflow_phase.set_value('temp', t_runoff, time)
        self.outflow_phase.set_value('flowrate', m_dot_2, time)
        self.outflow_phase.set_value('quality', quality, time)

        self.state_phase.add_row(time, turbine)

        self.state_phase.set_value('power', power, time)

        return time
