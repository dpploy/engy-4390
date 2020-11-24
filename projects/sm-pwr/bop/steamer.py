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

        self.primary_volume = 15*unit.meter**3
        self.secondary_volume = 25*unit.meter**3
        #self.ht_coeff = 4416194*unit.watt/unit.kelvin
        self.ht_coeff = 41720000.29*unit.watt/unit.kelvin

        # Initialization
        self.primary_mass_dens = 1*unit.gram/unit.cc
        self.primary_cp = 4.298 * unit.kj/unit.kg/unit.kelvin
        self.tau_primary = 1*unit.minute
        self.t_sat = 516 #K

        self.primary_inflow_pressure = 12.8*unit.mega*unit.pascal
        self.primary_inflow_temp = (30+273.15)*unit.kelvin
        self.primary_inflow_mass_flowrate = 666

        self.secondary_mass_dens = 1*unit.gram/unit.cc
        self.secondary_cp = 4.298 * unit.kj/unit.kg/unit.kelvin
        self.tau_secondary = 1*unit.minute

        self.secondary_inflow_pressure = 3.4*unit.mega*unit.pascal
        self.secondary_inflow_temp = (20+273.15)*unit.kelvin
        self.secondary_inflow_mass_flowrate = 67

        self.primary_outflow_temp = (20 + 273.15)*unit.kelvin
        self.primary_outflow_mass_flowrate = 67
        self.primary_outflow_pressure = 12.8*unit.mega*unit.pascal

        self.secondary_outflow_mass_flowrate = 67
        self.secondary_outflow_temp = (20 + 273.15)*unit.kelvin
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

            self.primary_inflow_temp = primary_inflow['temperature']
            self.primary_inflow_pressure = primary_inflow['pressure']
            self.primary_inflow_mass_flowrate = primary_inflow['mass_flowrate']

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
            secondary_outflow['mass_flowrate'] = flowrate

            self.send((msg_time, secondary_outflow), 'secondary-outflow')

    def __step(self, time=0.0):
        """ODE IVP problem.
        """

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

        temp_p = u_vec[0] # primary outflow temp
        temp_s = u_vec[1] # secondary outflow temp

        # Update phases
        primary_outflow = self.primary_outflow_phase.get_row(time)
        secondary_outflow = self.secondary_outflow_phase.get_row(time)

        time += self.time_step

        self.primary_outflow_phase.add_row(time, primary_outflow)
        self.primary_outflow_phase.set_value('temp', temp_p, time)

        self.secondary_outflow_phase.add_row(time, secondary_outflow)
        self.secondary_outflow_phase.set_value('temp', temp_s, time)
        self.secondary_outflow_phase.set_value('flowrate',
                                               self.secondary_inflow_mass_flowrate,
                                               time)
        self.secondary_outflow_phase.set_value('pressure',
                                               self.secondary_outflow_pressure,
                                               

                                               time)

        return time

    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.
           Neutron density, delayed neutron emmiter concentrations,
           termperature of fuel, and temperature of coolant.
        """

        u_vec = np.empty(0, dtype=np.float64)

        temp_p = self.primary_outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp_p)

        temp_s = self.secondary_outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp_s)

        return u_vec

    def __f_vec(self, u_vec, time):

        temp_p = u_vec[0] # get temperature of primary outflow

        temp_s = u_vec[1] # get temperature of secondary outflow

        # initialize f_vec to zero
        f_tmp = np.zeros(2, dtype=np.float64) # vector for f_vec return

        #-----------------------
        # primary energy balance
        #-----------------------
        rho_p = 1/(steam_table._Region1(self.primary_inflow_temp,self.primary_inflow_pressure/unit.mega)["v"])
        cp_p = steam_table._Region1(self.primary_inflow_temp,self.primary_inflow_pressure/unit.mega)["cp"]
        vol_p = self.primary_volume
        
        
        temp_p_in = self.primary_inflow_temp

        tau_p = self.tau_primary

        #-----------------------
        # secondary energy balance
        #-----------------------
        rho_s = 1/(steam_table._Region1((self.secondary_inflow_temp),self.secondary_inflow_pressure/unit.mega)["v"])
        cp_s = steam_table._Region1((self.secondary_inflow_temp),self.secondary_inflow_pressure/unit.mega)["cp"]
        vol_s = self.secondary_volume

        temp_s_in = self.secondary_inflow_temp

        tau_s = self.tau_secondary
        
        #-----------------------
        # calculations
        #-----------------------
        heat_sink = self.__heat_sink_rate(temp_p_in, temp_s_in, cp_p, cp_s)

        f_tmp[0] = - 1/tau_p * (temp_p - temp_p_in) + 1./rho_p/cp_p/vol_p * heat_sink
        
        
        


        heat_source = - heat_sink
        
        temp_s_out = - 1/tau_s * (temp_s - temp_s_in) + 1./rho_s/cp_s/vol_s * heat_source

        q_total = (temp_s_out-temp_s_in)*cp_s
        q_heat = (self.t_sat - temp_s_in)*cp_s
       
        if q_total < q_heat:
            f_tmp[1] = temp_s_out
            
        else:
            h_v = steam_table._Region4(self.secondary_inflow_pressure,1)['h']
            h_l = steam_table._Region4(self.secondary_inflow_pressure,0)['h']
            h_vap = h_v-h_l
            
            quality = (h_vap - (q_total - q_heat))/h_vap
            f_tmp[1] = steam_table._Region4(self.secondary_inflow_pressure,quality)['T']
        
        
        #print(temp_p_in,f_tmp[0], f_tmp[1])
        #f_tmp[1] = - 1/tau_s * (temp_s - temp_s_in) + 1./rho_s/cp_s/vol_s * heat_source

        return f_tmp

    def __heat_sink_rate(self, temp_p, temp_s, cp_p, cp_s):
        """Cooling rate of primary."""

        c_min = cp_s*self.secondary_inflow_mass_flowrate
        
        #ntu = self.ht_coeff/(c_min)
        ntu = 4.5
        
        c_r = (c_min)/(cp_p*self.primary_inflow_mass_flowrate)
        
        
        
        eta = (1-np.exp(-ntu*(1-c_r)))/(1-c_r*np.exp(-ntu*(1-c_r)))
        
       
        q_p = - eta*c_min*(temp_p - temp_s)
        
        


        return q_p
