#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import unit

import iapws.iapws97 as steam_table

import math
from scipy.integrate import odeint
import numpy as np

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
        self.tau = 1.5
        self.volume = 5

        # Initialization
        self.inflow_pressure = 0.008066866
        self.inflow_temp = 20+273
        self.inflow_mass_flowrate = 67
        self.inflow_quality = 0

        self.outflow_temp = 20 + 273.15
        self.outflow_mass_flowrate = 67
        self.outflow_pressure = 2
        
        self.heat_source_rate = -50000 #W

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
                         formal_name='P_2', unit='MPa',
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
        # one way "from" turbine

        # receive from
        if self.get_port('inflow').connected_port:

            self.send(time, 'inflow')

            (check_time, inflow) = self.recv('inflow')
            print(check_time,time)
            assert abs(check_time-time) <= 1e-6

            self.inflow_temp = inflow['temperature']
            self.inflow_pressure = inflow['pressure']
            self.inflow_mass_flowrate = inflow['mass_flowrate']

        # Interactions in the feed water inflow port
        #-----------------------------------------
        # one way "to" feedwater

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

        # Get state values
        t_exit = self.inflow_temp
        p_out = steam_table._PSat_T(self.inflow_temp)
        flow_out = self.inflow_mass_flowrate
        # Get state values
        u_0 = self.__get_state_vector(time)

        t_interval_sec = np.linspace(time, time+self.time_step, num=2)

        max_n_steps_per_time_step = 1000 # max number of nonlinear algebraic solver
                                         # iterations per time step

        (u_vec_hist, info_dict) = odeint(self.__f_vec, u_0, t_interval_sec,
                                         rtol=1e-4, atol=1e-8,
                                         mxstep=max_n_steps_per_time_step,
                                         full_output=True, tfirst=False)

        assert info_dict['message'] == 'Integration successful.', info_dict['message']

        u_vec = u_vec_hist[1, :]  # solution vector at final time step

        temp = u_vec[0] # primary outflow temp

        #update state variables
        condenser_outflow = self.outflow_phase.get_row(time)

        time += self.time_step

        self.outflow_phase.add_row(time, condenser_outflow)
        self.outflow_phase.set_value('temp', t_exit, time)
        self.outflow_phase.set_value('flowrate', flow_out, time)
        self.outflow_phase.set_value('pressure', p_out, time)

        return time
    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.

        """

        u_vec = np.empty(0, dtype=np.float64)

        temp = self.outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp)

        return  u_vec

    def __f_vec(self, u_vec, time):

        temp = u_vec[0]

        # initialize f_vec to zero
        f_tmp = np.zeros(1, dtype=np.float64) # vector for f_vec return

        #-----------------------
        # primary energy balance
        #-----------------------
        if self.inflow_quality == 0:
            rho = 1/(steam_table._Region1(self.inflow_temp,self.inflow_pressure)["v"])
            cp = steam_table._Region1(self.inflow_temp,self.inflow_pressure)["cp"]
        elif 0 < self.inflow_quality < 1:
            rho = 1/(steam_table._Region4(self.inflow_pressure,self.inflow_quality)["v"])
            cp = steam_table._Region4(self.inflow_pressure, self.inflow_quality)["cp"]
        else:
            rho = 1/(steam_table._Region2(self.inflow_temp,self.inflow_pressure)["v"])
            cp = steam_table._Region2(self.inflow_temp,self.inflow_pressure)["cp"]
            
        
        vol = self.volume

        temp_in = self.inflow_temp

        tau = self.tau

        #-----------------------
        # calculations
        #-----------------------
        heat_source = self.heat_source_rate

        f_tmp[0] = - 1/tau * (temp - temp_in) + 1./rho/cp/vol * heat_source

        return f_tmp