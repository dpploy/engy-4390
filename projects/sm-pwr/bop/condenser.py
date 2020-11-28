#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import unit

from iapws import IAPWS97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class Condenser(Module):
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

        self.port_names_expected = ['inflow', 'outflow']

        # Units

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

        self.inflow_pressure = 34*unit.bar
        self.inflow_temp = (20+273)*unit.K
        self.inflow_mass_flowrate = 67*unit.kg/unit.second

        self.outflow_temp = (20+273.15)*unit.K
        self.outflow_mass_flowrate = 67*unit.kg/unit.second
        self.outflow_pressure = 34*unit.bar

        # Outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Condenser Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T_2$',
                        info='Condenser Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='Pa',
                         value=self.outflow_pressure,
                         latex_name=r'$P_2$',
                         info='Condenser Outflow Pressure')

        quantities.append(press)

        self.outflow_phase = Phase(time_stamp=self.initial_time,
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

        # Interactions in the outflow port
        #-----------------------------------------
        # one way "to" outflow

        # send to
        if self.get_port('outflow').connected_port:

            msg_time = self.recv('outflow')
            assert msg_time <= time

            temp = self.outflow_phase.get_value('temp', msg_time)
            pressure = self.outflow_phase.get_value('pressure', msg_time)
            outflow = dict()
            outflow['temperature'] = temp
            outflow['pressure'] = pressure
            outflow['mass_flowrate'] = self.outflow_mass_flowrate
            self.send((msg_time, outflow), 'outflow')

    def __step(self, time=0.0):

        # Comput state
        #print('inflow pressure [bar] = ',self.inflow_pressure/unit.bar)
        water = steam_table(P=self.inflow_pressure/unit.mega/unit.pascal, x=0.0)

        outflow_temp = water.T
        outflow_mass_flowrate = self.inflow_mass_flowrate
        outflow_pressure = self.inflow_pressure

        #update state variables
        condenser_outflow = self.outflow_phase.get_row(time)

        time += self.time_step

        self.outflow_phase.add_row(time, condenser_outflow)
        self.outflow_phase.set_value('temp', outflow_temp, time)
        self.outflow_phase.set_value('flowrate', outflow_mass_flowrate , time)
        self.outflow_phase.set_value('pressure', outflow_pressure, time)

        return time
