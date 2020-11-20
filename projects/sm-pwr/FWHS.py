#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import scipy.constants as unit

import iapws.iapws97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class FWHS(Module):
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

        self.port_names_expected = ['inturb','incond', 'outflow']

        # Units
        unit.kg = unit.kilo*unit.gram
        unit.meter = 1.0
        unit.cm = unit.centi*unit.meter
        unit.second = 1.0
        unit.pascal = 1.0
        unit.joule = 1.0
        unit.kj = unit.kilo*unit.joule
        unit.kelvin = 1.0
        unit.watt = 1.0

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

        self.inturb_pressure = 0.0
        self.inturb_temp = 0.0
        self.inturb_mass_flowrate = 0.0
        
        self.incond_pressure = 0.0
        self.incond_temp = 0.0
        self.incond_mass_flowrate = 0.0

        self.outflow_temp = 20 + 273.15
        self.outflow_temp_ss = 422 #k
        self.outflow_mass_flowrate = 0.0
        self.outflow_pressure = 0.0
        self.outflow_pressure_ss = 3.44738 #MPa
        self.outflow_h = steam_table._Region4(p_out, 0)
        
        # Outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T_2$',
                        info='Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='MPa',
                         value=self.outflow_pressure,
                         latex_name=r'$P_2$',
                         info='Outflow Pressure')

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

        # Interactions in the feed water inflow port
        #-----------------------------------------
        # one way "to" steam geneator

        # send to
        if self.get_port('outflow').connected_port:

            msg_time = self.recv('outflow')
            assert msg_time <= time

            temp = self.outflow_phase.get_value('temp', msg_time)
            pressure = self.outflow_phase.get_value('pressure', msg_time)
            outflow = dict()
            outflow['temperature'] = temp
            outflow['pressure'] = pressure
            outflow['flowrate'] = self.outflow_mass_flowrate
            self.send((msg_time, outflow), 'outflow')

        # Interactions in the incond port
        #----------------------------------------
        # one way "from" condenser

        # receive from
        if self.get_port('incond').connected_port:

            self.send(time, 'incond')

            (check_time, incond) = self.recv('incond')
            assert abs(check_time-time) <= 1e-6

            self.incond_temp = incond['temperature']
            self.incond_pressure = incond['pressure']
            self.incond_mass_flowrate = incond['mass_flowrate']
            
         # Interactions in the inturb port
        #----------------------------------------
        # one way "from" turbine

        # receive from
        if self.get_port('inturb').connected_port:

            self.send(time, 'inturb')

            (check_time, inturb) = self.recv('inturb')
            assert abs(check_time-time) <= 1e-6

            self.inturb_temp = inturb['temperature']
            self.inturb_pressure = inturb['pressure']
            self.inturb_mass_flowrate = inturb['mass_flowrate']
            
       

    def __step(self, time=0.0):

        # Get state values
        flowturb = self.inturb_mass_flowrate
        flowcond = self.incond_mass_flowrate
        hturb_in = steam_table._Region1(self.inturb_temp,self.inturb_pressure)['h']
        hcond_in = steam_table._Region1(self.incond_temp,self.incond_pressure)['h']
        t_exit = self.outflow_temp_ss
        p_out = self.outflow_pressure_ss
        flow_out = flowturb+flowcond
        h_exit = steam_table._Region1(t_out,p_out)['h']
        q_removed = flow_out*h_exit-flowturb*hturbin-flowcond*hcondin

        #update state variables
        condenser_outflow = self.outflow_phase.get_row(time)

        time += self.time_step

        self.outflow_phase.add_row(time, flow_out)
        self.outflow_phase.set_value('temp', t_exit, time)
        self.outflow_phase.set_value('flowrate', flow_out, time)
        self.outflow_phase.set_value('pressure', p_out, time)

        return time
