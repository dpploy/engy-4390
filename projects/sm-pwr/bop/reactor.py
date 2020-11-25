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

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class SMPWR(Module):
    """Small modular pressurized boiling water single-point reactor.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: steam generator.
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['coolant-inflow', 'coolant-outflow']

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters

        # Data pertaining to one-group energy neutron balance
        self.gen_time = 1.0e-4*unit.second
        self.beta = 5.9e-3
        self.diff_coeff = 0.2939*unit.cm
        self.k_infty = 1.49826
        self.buckling = 1.538e-4

        self.alpha_n = -3e-6 # control rod reactivity worth

        self.n_dens_ss_operation = 5e14/2200/unit.meter**3

        #Delayed neutron emission
        self.species_decay = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01] # 1/sec
        self.species_rel_yield = [0.033, 0.219, 0.196, 0.395, 0.115, 0.042]

        #Data pertaining to two-temperature heat balances
        self.fis_energy = 180 * 1.602e-13 *unit.joule # J/fission
        self.sigma_f_o = 586.2 * 100 * 1e-30 * unit.meter**2
        self.temp_o = 20 + 273.15 # K
        self.temp_c_ss_operation = 265 + 273.15# K desired ss operation temp of coolant
        self.thermal_neutron_velo = 2200*unit.meter/unit.second

        self.fis_nuclide_num_dens = 9.84e26/unit.meter**3 # (fissile nuclei)/m3
        self.condenser_pressure = 0.008066866 #MPa
        self.temp_inlet_ss = 265 + 273.15

        self.fuel_dens = 10800*unit.kg/unit.meter**3
        self.cp_fuel = 300*unit.joule/unit.kg/unit.kelvin
        self.fuel_volume = .8565*unit.meter**3

        self.coolant_mass_flowrate = 666 #kg/s
        self.coolant_dens = 669.2294308156266*unit.kg/unit.meter**3
        self.cp_coolant = 1000*5.382268683703659 # J/(mol K) - > J/(kg K)
        self.coolant_volume = 2.8*unit.meter**3
        self.coolant_pressure = 12.8 #MPa

        self.ht_coeff = 1300000*unit.watt/unit.kelvin

        self.tau = 2.8*unit.second # coolant flow residence time

        # Initialization
        self.n_dens_ref = 1.0
        self.q_0 = 1./self.gen_time # pulse neutron source
        rho_0_over_beta = 0.35# $

        self.n_0 = 0.0 # neutronless steady state before start up
        self.rho_0 = rho_0_over_beta * self.beta # neutron source

        lambda_vec = np.array(self.species_decay, dtype=np.float64)
        beta_vec = np.array(self.species_rel_yield, dtype=np.float64) * self.beta
        c_vec_0 = beta_vec/lambda_vec/self.gen_time * self.n_0

        self.temp_f_0 = self.temp_o
        self.temp_c_0 = self.temp_o

        self.inflow_cool_temp = self.temp_o

        # Coolant outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_c', unit='kg/s',
                            value=self.coolant_mass_flowrate,
                            latex_name=r'$q_c$',
                            info='Reactor Outflow Coolant Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_c', unit='K',
                        value=self.temp_c_0,
                        latex_name=r'$T_c$',
                        info='Reactor Outflow Coolant Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_c', unit='Pa',
                         latex_name=r'$P_c$',
                         info='Reactor Outflow Coolant Pressure')

        quantities.append(press)

        self.coolant_outflow_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s',
                                           quantities=quantities)

        # Neutron phase history
        quantities = list()

        neutron_dens = Quantity(name='neutron-dens',
                                formal_name='n', unit='1/m^3',
                                value=self.n_0,
                                latex_name=r'$n$',
                                info='Reactor Neutron Density')

        quantities.append(neutron_dens)

        delayed_neutron_cc = Quantity(name='delayed-neutrons-cc',
                                      formal_name='c_i', unit='1/m^3 ',
                                      value=c_vec_0,
                                      latex_name=r'$c_i$',
                                      info='Reactor Delayed Neutron Precursors')

        quantities.append(delayed_neutron_cc)

        self.neutron_phase = Phase(time_stamp=self.initial_time, time_unit='s',
                                   quantities=quantities)

        # Reactor phase
        quantities = list()

        fuel_temp = Quantity(name='fuel-temp',
                             formal_name='T_f', unit='K',
                             value=self.temp_f_0,
                             latex_name=r'$T_f$',
                             info='Nuclear Fuel Temperature')
        quantities.append(fuel_temp)
        
        inlet_temp = Quantity(name='inlet-temp',
                             formal_name='T_in', unit='K',
                             value=self.inflow_cool_temp,
                             latex_name=r'$T_{in}$',
                             info='Inflow Coolant Temperature')
        quantities.append(inlet_temp)

        self.reactor_phase = Phase(time_stamp=self.initial_time, time_unit='s',
                                   quantities=quantities)


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

        # Interactions in the coolant-outflow port
        #-----------------------------------------
        # One way "to" coolant-outflow

        # Send to
        if self.get_port('coolant-outflow').connected_port:

            msg_time = self.recv('coolant-outflow')
            assert msg_time <= time

            outflow_cool_temp = self.coolant_outflow_phase.get_value('temp', msg_time)
            coolant_outflow = dict()
            coolant_outflow['temperature'] = outflow_cool_temp
            coolant_outflow['pressure'] = self.coolant_pressure
            coolant_outflow['mass_flowrate'] = self.coolant_mass_flowrate
            self.send((msg_time, coolant_outflow), 'coolant-outflow')

        # Interactions in the coolant-inflow port
        #----------------------------------------
        # One way "from" coolant-inflow

        # Receive from
        if self.get_port('coolant-inflow').connected_port:

            self.send(time, 'coolant-inflow')

            (check_time, inflow_coolant) = self.recv('coolant-inflow')
            assert abs(check_time-time) <= 1e-6

            self.inflow_cool_temp = inflow_coolant['temperature']
            self.coolant_mass_flowrate = inflow_coolant['mass_flowrate']
            self.coolant_pressure = inflow_coolant['pressure']

    def __step(self, time=0.0):
        r"""ODE IVP problem.
        Given the initial data at :math:`t=0`,
        :math:`u = (u_1(0),u_2(0),\ldots)`
        solve :math:`\frac{\text{d}u}{\text{d}t} = f(u)` in the interval
        :math:`0\le t \le t_f`.

        Parameters
        ----------
        time: float
            Time in SI unit.

        Returns
        -------
        None
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

        n_dens = u_vec[0]
        c_vec = u_vec[1:7]
        fuel_temp = u_vec[7]
        cool_temp = u_vec[8]

        # Update phases
        coolant_outflow = self.coolant_outflow_phase.get_row(time)
        neutrons = self.neutron_phase.get_row(time)
        reactor = self.reactor_phase.get_row(time)

        time += self.time_step

        self.coolant_outflow_phase.add_row(time, coolant_outflow)
        self.neutron_phase.add_row(time, neutrons)
        self.reactor_phase.add_row(time, reactor)

        self.coolant_outflow_phase.set_value('temp', cool_temp, time)
        self.neutron_phase.set_value('neutron-dens', n_dens, time)
        self.neutron_phase.set_value('delayed-neutrons-cc', c_vec, time)
        self.reactor_phase.set_value('fuel-temp', fuel_temp, time)
        self.reactor_phase.set_value('inlet-temp', self.inflow_cool_temp, time)
        return time

    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.
           Neutron density, delayed neutron emmiter concentrations,
           termperature of fuel, and temperature of coolant.
        """

        u_vec = np.empty(0, dtype=np.float64)

        neutron_dens = self.neutron_phase.get_value('neutron-dens', time)
        u_vec = np.append(u_vec, neutron_dens)

        delayed_neutrons_cc = self.neutron_phase.get_value('delayed-neutrons-cc', time)
        u_vec = np.append(u_vec, delayed_neutrons_cc)

        fuel_temp = self.reactor_phase.get_value('fuel-temp', time)
        u_vec = np.append(u_vec, fuel_temp)

        temp = self.coolant_outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp)

        return u_vec

    def __alpha_tn_func(self, temp):
        """Single energy group formula.
        """

        B2 = self.buckling
        D = self.diff_coeff
        k_infty = self.k_infty
        Ea = .022  #/cm
        To = self.temp_o

        alpha_tn = -1.0 / 2.0 * B2 * D / (k_infty * Ea * math.sqrt(To * temp))

        return alpha_tn

    def __rho_func(self, time, n_dens, temp):
        """Reactivity function.

        Parameters
        ----------
        t: float
            Time.
        temp_f: float
            Temperature at time t.
        params: dict
            Dictionary of quantities.

        Returns
        -------
        rho_t: float
            Value of reactivity.

        Examples
        --------
        """

        temp_ref = self.temp_c_ss_operation
        alpha_n = self.alpha_n
        n_dens_ref = self.n_dens_ref

        alpha_tn = self.__alpha_tn_func(temp)

        rho_t = self.rho_0 + alpha_n * (n_dens - n_dens_ref) + \
                             alpha_tn * (temp - temp_ref)

        return rho_t

    def __q_source(self, time):
        """Neutron source delta function.

        Parameters
        ----------
        t: float
            Time.

        Returns
        -------
        q: float
            Value of source.

        Examples
        --------
        """

        #broken on Mac/Windows if time <= 100*unit.milli*unit.second: # small time value
        if time <= 100*unit.milli: # small time value
            qval = self.q_0
        else:
            qval = 0.0

        return qval

    def __sigma_fis_func(self, temp):
        """Effective microscopic fission cross section.
        """

        sigma_f = self.sigma_f_o * math.sqrt(self.temp_o/temp) * \
                  math.sqrt(math.pi)/2.0

        return sigma_f

    def __nuclear_pwr_dens_func(self, time, temp, n_dens):
        """Place holder for implementation.
        """

        rxn_heat = self.fis_energy # get fission reaction energy J per reaction

        sigma_f = self.__sigma_fis_func(temp)

        fis_nuclide_num_dens = self.fis_nuclide_num_dens #  #/m3

        macro_sigma_f = sigma_f * fis_nuclide_num_dens # macroscopic cross section

        v_o = self.thermal_neutron_velo # m/s

        neutron_flux = n_dens * self.n_dens_ss_operation * v_o

        #reaction rate density
        rxn_rate_dens = macro_sigma_f * neutron_flux

        # nuclear heating power density
        q3prime = - rxn_heat * rxn_rate_dens # exothermic reaction W/m3)

        return q3prime

    def __heat_sink_rate(self, temp_f, temp_c):
        """Cooling rate."""

        ht_coeff = self.ht_coeff

        q_f = - ht_coeff * (temp_f - temp_c)

        return q_f

    def __f_vec(self, u_vec, time):

        n_dens = u_vec[0] # get neutron dens

        c_vec = u_vec[1:-2] # get delayed neutron emitter concentration

        temp_f = u_vec[-2] # get temperature of fuel

        temp_c = u_vec[-1] # get temperature of coolant

        # initialize f_vec to zero
        lambda_vec = np.array(self.species_decay, dtype=np.float64)
        n_species = len(lambda_vec)

        f_tmp = np.zeros(1+n_species+2, dtype=np.float64) # vector for f_vec return

        #----------------
        # neutron balance
        #----------------

        rho_t = self.__rho_func(time, n_dens, (temp_f+temp_c)/2.0)

        beta = self.beta
        gen_time = self.gen_time

        assert len(lambda_vec) == len(c_vec)

        q_source_t = self.__q_source(time)

        f_tmp[0] = (rho_t - beta)/gen_time * n_dens + lambda_vec @ c_vec + q_source_t

        #-----------------------------------
        # n species balances (implicit loop)
        #-----------------------------------
        species_rel_yield = self.species_rel_yield
        beta_vec = np.array(species_rel_yield) * beta

        assert len(beta_vec) == len(c_vec)

        f_tmp[1:-2] = beta_vec / gen_time * n_dens - lambda_vec * c_vec

        #--------------------
        # fuel energy balance
        #--------------------
        rho_f = self.fuel_dens
        cp_f = self.cp_fuel
        vol_fuel = self.fuel_volume

        pwr_dens = self.__nuclear_pwr_dens_func(time, (temp_f+temp_c)/2, n_dens)

        heat_sink = self.__heat_sink_rate(temp_f, temp_c)

        f_tmp[-2] = -1/rho_f/cp_f * (pwr_dens - heat_sink/vol_fuel)

        #-----------------------
        # coolant energy balance
        #-----------------------
        rho_c = self.coolant_dens
        cp_c = self.cp_coolant
        vol_cool = self.coolant_volume

        temp_in = self.inflow_cool_temp

        tau = self.tau

        heat_source = - heat_sink

        f_tmp[-1] = - 1/tau * (temp_c - temp_in) + 1./rho_c/cp_c/vol_cool * heat_source

        return f_tmp
