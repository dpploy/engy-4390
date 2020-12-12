#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module.
   PWR for NuScale BOP.

       + Coolant outflow temperature: 543 F (283.9 C) (w/ a 100 F rise)
       + Coolant inflow temperature: 497 F (258.3 C)
       + Coolant pressure: 1850 psi (127.6 bar)
       + Coolant mass flowrate:
             - Avg: 4.66e+6 (avg) lb/h
             - Max: 5.24e6
             - Min: 4.27e6
       + Coolant mass flux: 0.49e+6 lb/h-ft^2
       + Fuel temperature: 1200 C
       + Core average centerline pellet temperature (2.5 kw/ft): 1375 F (746.11 C)
       + Core rod centerline pellet temperature (6.5 kw/ft): 2075 F (1135 C)
       + Heat transfer area on fuel: 6275 ft^2
       + Core flow area: 9.79 ft^2 (0.91 m^2)
       + Core flow diameter: 1.08 m
       + Core fuel volume: 0.841 m^3
       + Core fuel mass: 9221.88 kg (UO2)
       + Core fuel mass density at 273.15K:  9221.88/0.841
       + Effective fuel length = 95.89 in (2.44 m)
       + Equivalent diameter of active core = 59.28 in
       + Power: 160 MW
       + Normal operation peak heat flux = 0.171e+6 Btu/hr-ft^2
       + Normal operation core average heat flux = 84044 Btu/hr-ft^2 (0.26512442 MW/m^2)
"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

import unit

from iapws import IAPWS97 as WaterProps

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
        self.save = True

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        self.shutdown = (False, 10*unit.minute)

        # Configuration parameters
        self.discard_tau_recording_before = 2*unit.minute

        # Data pertaining to one-group energy neutron balance
        self.gen_time = 1.0e-4*unit.second
        self.beta = 5.9e-3
        self.diff_coeff = 0.2939*unit.cm
        self.k_infty = 1.49826
        self.buckling = 1.538e-4

        #self.alpha_n = -3e-6 # control rod reactivity worth
        self.alpha_n = -1e-5 # control rod reactivity worth

        self.n_dens_ss_operation = 2*5e14/2200/unit.meter**3

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

        self.core_dens = 10800*unit.kg/unit.meter**3
        self.cp_core = 300*unit.joule/unit.kg/unit.kelvin
        self.core_volume = 0.841*unit.meter**3

        self.coolant_mass_flowrate_ss = 4.66e6*unit.lb/unit.hour
        self.flowrate_relaxation_startup = 5*unit.minute
        self.flowrate_relaxation_shutdown = 1*unit.minute

        self.coolant_pressure = 1850*unit.psi  # altered by coolant-inflow port

        self.core_heat_transfer_area = 6275*unit.foot**2
        self.core_flow_area = 9.79*unit.foot**2
        self.core_length = 95.89*unit.inch
        self.cladding_k = 12.1*unit.watt/unit.meter/unit.kelvin
        self.coolant_volume = self.core_flow_area * self.core_length

        self.coolant_quality = 0.0 # not used; to be used in heat flux correlations

        # Initialization
        self.n_dens_ref = 1.0
        self.q_0 = 1./self.gen_time # pulse neutron source
        rho_0_over_beta = 0.25# $

        self.n_0 = 0.0 # neutronless steady state before start up
        self.rho_0 = rho_0_over_beta * self.beta # neutron source

        lambda_vec = np.array(self.species_decay, dtype=np.float64)
        beta_vec = np.array(self.species_rel_yield, dtype=np.float64) * self.beta
        c_vec_0 = beta_vec/lambda_vec/self.gen_time * self.n_0

        self.temp_f_0 = self.temp_o + 20*unit.K
        self.temp_c_0 = self.temp_o + 20*unit.K

        self.inflow_cool_temp = self.temp_o + 20*unit.K

        # Derived quantities

        # Coolant outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='m_c', unit='kg/s',
                            value=0.0,
                            latex_name=r'$\dot{m}_c$',
                            info='Reactor Coolant Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_c', unit='K',
                        value=self.temp_c_0,
                        latex_name=r'$T_c$',
                        info='Reactor Coolant Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_c', unit='Pa',
                         value=self.coolant_pressure,
                         latex_name=r'$P_c$',
                         info='Reactor Coolant Outflow Pressure')

        quantities.append(press)

        press = Quantity(name='quality',
                         formal_name='X_c', unit='%',
                         value=0.0,
                         latex_name=r'$\chi_c$',
                         info='Reactor Coolant Outflow Quality')

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

        # State phase
        quantities = list()

        core_temp = Quantity(name='core-temp',
                             formal_name='T_f', unit='K',
                             value=self.temp_f_0,
                             latex_name=r'$T_f$',
                             info='Reactor Core Temperature')
        quantities.append(core_temp)

        inlet_temp = Quantity(name='inlet-temp',
                              formal_name='T_in', unit='K',
                              value=self.inflow_cool_temp,
                              latex_name=r'$T_{in}$',
                              info='Reactor Coolant Inflow Temperature')
        quantities.append(inlet_temp)

        pwr = Quantity(name='power',
                       formal_name='Pth', unit='W',
                       value=0.0,
                       latex_name=r'$P_{th}$',
                       info='Reactor Power')
        quantities.append(pwr)

        rey = Quantity(name='reynolds',
                       formal_name='Rey_c', unit='',
                       value=0.0,
                       latex_name=r'$R_{c}$',
                       info='Reactor Coolant Reynolds Number')
        quantities.append(rey)

        water = WaterProps(T=self.temp_c_0, P=self.coolant_pressure/unit.mega/unit.pascal)

        prtl = Quantity(name='prandtl',
                        formal_name='Pr_c', unit='',
                        value=water.Prandt,
                        latex_name=r'$Pr_{c}$',
                        info='Reactor Coolant Prandtl Number')
        quantities.append(prtl)

        q2prime = Quantity(name='heatflux',
                           formal_name="q''", unit=r'W/m$^2$',
                           value=0.0,
                           latex_name=r"$q''$",
                           info='Reactor Core Average Heat Flux')
        quantities.append(q2prime)

        nusselt = Quantity(name='nusselt',
                           formal_name='Nu_c', unit='',
                           value=0.0,
                           latex_name=r'$Nu_{c}$',
                           info='Reactor Coolant Nusselt Number')
        quantities.append(nusselt)

        tau = Quantity(name='tau',
                       formal_name='tau', unit='s',
                       value=0.0,
                       latex_name=r'$\tau_{c}$',
                       info='Reactor Coolant Residence Time')
        quantities.append(tau)

        self.state_phase = Phase(time_stamp=self.initial_time, time_unit='s',
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

            # Evolve one time step
            #---------------------
            time = self.__step(time)

            # Communicate information
            #------------------------
            self.__call_ports(time)

            # Shudown procedure
            #------------------
            if self.shutdown[0] and time > self.shutdown[1]:
                self.rho_0 = -1.0
                self.species_rel_yield = len(self.species_rel_yield)*[0.0]

    def __call_ports(self, time):

        # Interactions in the coolant-outflow port
        #-----------------------------------------
        # One way "to" coolant-outflow

        # Send to
        if self.get_port('coolant-outflow').connected_port:

            msg_time = self.recv('coolant-outflow')
            assert msg_time <= time

            temp = self.coolant_outflow_phase.get_value('temp', msg_time)
            press = self.coolant_outflow_phase.get_value('pressure', msg_time)
            flowrate = self.coolant_outflow_phase.get_value('flowrate', msg_time)
            chi = self.coolant_outflow_phase.get_value('quality', msg_time)

            coolant_outflow = dict()
            coolant_outflow['temperature'] = temp
            coolant_outflow['pressure'] = press
            coolant_outflow['mass_flowrate'] = flowrate
            coolant_outflow['quality'] = chi

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
            self.coolant_quality = inflow_coolant['quality']

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

        max_n_steps_per_time_step = 1500 # max number of nonlinear algebraic solver
                                         # iterations per time step

        (u_vec_hist, info_dict) = odeint(self.__f_vec, u_0, t_interval_sec,
                                         rtol=1e-7, atol=1e-8,
                                         mxstep=max_n_steps_per_time_step,
                                         full_output=True, tfirst=False)

        assert info_dict['message'] == 'Integration successful.', info_dict['message']

        u_vec = u_vec_hist[1, :]  # solution vector at final time step

        n_dens = u_vec[0]
        c_vec = u_vec[1:7]
        mass_flowrate = u_vec[7]
        core_temp = u_vec[8]
        cool_temp = u_vec[9]

        # Update phases
        coolant_outflow = self.coolant_outflow_phase.get_row(time)
        neutrons = self.neutron_phase.get_row(time)
        reactor = self.state_phase.get_row(time)

        time += self.time_step

        self.coolant_outflow_phase.add_row(time, coolant_outflow)
        self.neutron_phase.add_row(time, neutrons)
        self.state_phase.add_row(time, reactor)

        self.coolant_outflow_phase.set_value('temp', cool_temp, time)
        self.coolant_outflow_phase.set_value('flowrate', mass_flowrate, time)

        self.neutron_phase.set_value('neutron-dens', n_dens, time)
        self.neutron_phase.set_value('delayed-neutrons-cc', c_vec, time)

        self.state_phase.set_value('core-temp', core_temp, time)
        self.state_phase.set_value('inlet-temp', self.inflow_cool_temp, time)

        # Coolant properties
        water = WaterProps(T=cool_temp, P=self.coolant_pressure/unit.mega/unit.pascal)
        if water.phase == 'Two phases':
            qual = water.x
            assert qual <= 0.4 # limit to low quality
            #rho_c = (1-qual)*water.Liquid.rho + qual*water.Vapor.rho
            rho_c = water.rho
            cp_c = (1-qual)*water.Liquid.cp + qual*water.Vapor.cp
            cp_c *= unit.kj/unit.kg/unit.K
            mu_c = (1-qual)*water.Liquid.mu + qual*water.Vapor.mu
            k_c = (1-qual)*water.Liquid.k + qual*water.Vapor.k
            prtl_c = (1-qual)*water.Liquid.Prandt + qual*water.Vapor.Prandt
        elif water.phase == 'Liquid':
            rho_c = water.Liquid.rho
            cp_c = water.Liquid.cp*unit.kj/unit.kg/unit.K
            k_c = water.Liquid.k
            mu_c = water.Liquid.mu
            prtl_c = water.Liquid.Prandt
        else:
            assert False,'Vapor not allowed.'

        # Reactor power
        pwr = mass_flowrate*cp_c*(cool_temp-self.inflow_cool_temp)

        if pwr <= 0: # case when reactor is heated by the coolant inflow
            pwr = 0.0
        self.state_phase.set_value('power', abs(pwr), time)

        # Reynolds number
        diameter = (4*self.core_flow_area/math.pi)**.5
        rey_c = 4*mass_flowrate / mu_c / math.pi / diameter
        self.state_phase.set_value('reynolds', rey_c, time)

        # Prandtl number
        self.state_phase.set_value('prandtl', prtl_c, time)

        # Heat flux and Nusselt number
        (heat_sink_rate, nusselt) = self.__heat_sink_rate(core_temp, water,
                                                          mass_flowrate)

        q2prime = - heat_sink_rate/self.core_heat_transfer_area
        if q2prime < 0.0: # case when reactor is heated by the coolant inflow
            q2prime = 0.0
        self.state_phase.set_value('heatflux', q2prime, time)
        self.state_phase.set_value('nusselt', nusselt, time)

        heat_rate_transfered = - heat_sink_rate

        # Coolant quality 
        # Coolant temperature is likely below saturation but there is quality in
        # view of local nucleate boiling.
        if water.phase == 'Liquid':
            water_sat_l = WaterProps(P=self.coolant_pressure/unit.mega/unit.pascal, x=0)
            water_sat_v = WaterProps(P=self.coolant_pressure/unit.mega/unit.pascal, x=1)
            spfc_h_sat_l = water_sat_l.Liquid.h * unit.kj/unit.kg
            spfc_h_sat_v = water_sat_v.Vapor.h * unit.kj/unit.kg

            heat_rate_latent = (spfc_h_sat_v-spfc_h_sat_l)*mass_flowrate
            heat_rate_sensible = (water_sat_l.T-self.inflow_cool_temp)*cp_c*mass_flowrate
        #print(heat_rate_transfered/unit.mega, heat_rate_sensible/unit.mega, heat_rate_latent/unit.mega)
        #print((heat_rate_transfered-heat_rate_sensible)/(heat_rate_latent-heat_rate_sensible)*100)

            quality = (heat_rate_transfered-heat_rate_sensible)/\
                      (heat_rate_latent-heat_rate_sensible)*100

            quality = quality if 0<=quality<=100 else 0

        elif water.phase == 'Two phases':
            quality = water.x

        self.coolant_outflow_phase.set_value('quality', quality, time)

        # Coolant residence time
        q_vol = mass_flowrate/rho_c

        tau = self.coolant_volume/q_vol \
              if q_vol > 0 and time > self.discard_tau_recording_before else 0

        self.state_phase.set_value('tau', tau, time)

        return time

    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.
           Neutron density, delayed neutron emmiter concentrations,
           temperature of core, and temperature of coolant.
        """

        u_vec = np.empty(0, dtype=np.float64)

        neutron_dens = self.neutron_phase.get_value('neutron-dens', time)
        u_vec = np.append(u_vec, neutron_dens)

        delayed_neutrons_cc = self.neutron_phase.get_value('delayed-neutrons-cc', time)
        u_vec = np.append(u_vec, delayed_neutrons_cc)

        flowrate = self.coolant_outflow_phase.get_value('flowrate', time)
        u_vec = np.append(u_vec, flowrate)

        core_temp = self.state_phase.get_value('core-temp', time)
        u_vec = np.append(u_vec, core_temp)

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
        """Scaled nuclear power density.
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

    def __f_vec(self, u_vec, time):

        n_dens = u_vec[0] # get neutron dens

        c_vec = u_vec[1:-3] # get delayed neutron emitter concentration

        mass_flowrate = u_vec[-3]

        temp_f = u_vec[-2] # get temperature of core

        temp_c = u_vec[-1] # get temperature of coolant


        # initialize f_vec to zero
        lambda_vec = np.array(self.species_decay, dtype=np.float64)
        n_species = len(lambda_vec)

        f_tmp = np.zeros(1+n_species+3, dtype=np.float64) # vector for f_vec return

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

        f_tmp[1:-3] = beta_vec / gen_time * n_dens - lambda_vec * c_vec

        #----------------------------
        # mass flowrate ("buoyancy")
        #----------------------------

        if self.shutdown[0] and time > self.shutdown[1]:
            tau = self.flowrate_relaxation_shutdown
            self.coolant_mass_flowrate_ss = 1.0*unit.kg/unit.second
        else:
            tau = self.flowrate_relaxation_startup

        f_tmp[-3] = - 1/tau * (mass_flowrate - self.coolant_mass_flowrate_ss)

        #---------------------------
        # Heating power calculations
        #---------------------------
        water = WaterProps(T=temp_c, P=self.coolant_pressure/unit.mega/unit.pascal)
        assert water.phase != 'Vapour'

        # Compute the heating sink power
        (heat_sink_pwr, _) = self.__heat_sink_rate(temp_f, water, mass_flowrate)

        #--------------------
        # core energy balance
        #--------------------
        rho_f = self.core_dens
        cp_f = self.cp_core
        vol_core = self.core_volume

        nuclear_pwr_dens = self.__nuclear_pwr_dens_func(time, (temp_f+temp_c)/2, n_dens)

        heat_sink_pwr_dens = heat_sink_pwr/vol_core

        f_tmp[-2] = -1/rho_f/cp_f * (nuclear_pwr_dens - heat_sink_pwr_dens)

        #-----------------------
        # coolant energy balance
        #-----------------------
        if water.phase == 'Two phases':
            qual = water.x
            assert qual <= 0.4 # limit to low quality
            cp_c = (1-qual)*water.Liquid.cp + qual*water.Vapor.cp
            cp_c *= unit.kj/unit.kg/unit.K
            #rho_c = (1-qual)*water.Liquid.rho + qual*water.Vapor.rho
            rho_c = water.rho
        elif water.phase == 'Liquid':
            cp_c = water.Liquid.cp*unit.kj/unit.kg/unit.K
            rho_c = water.Liquid.rho
        else:
            assert False,'Vapor not allowed.'

        vol_cool = self.coolant_volume

        temp_in = self.inflow_cool_temp

        q_vol = mass_flowrate/rho_c
        if q_vol > 0:
            tau = vol_cool / q_vol
        else:
            tau = 1*unit.hour

        # Heating source power
        heat_source_pwr = - heat_sink_pwr
        heat_source_pwr_dens = heat_source_pwr/vol_cool

        f_tmp[-1] = - 1/tau * (temp_c - temp_in) + 1./rho_c/cp_c * heat_source_pwr_dens

        return f_tmp

    def __heat_sink_rate(self, temp_f, water, coolant_mass_flowrate):
        """Cooling rate of the core.

           Assumptions
           -----------

           + Coolant: overall ranging from one phase heat tranfer to transition
             nucleate boiling.
        """

        assert water.phase != 'Vapour'

        # Primary props
        temp_c = water.T

        #print(water.P*unit.mega/unit.bar)
        #print(temp_c-273.15)

        water_sat = WaterProps(P=water.P, x=0.0)
        temp_c_sat = water_sat.T

        # Overall condition on colant; locally there may be nucleate boiling
        assert temp_c <= temp_c_sat

        if water.phase == 'Two phases':
            qual = water.x
            assert qual <= 0.4 # limit to low quality
            cp_c = (1-qual)*water.Liquid.cp + qual*water.Vapor.cp
            cp_c *= unit.kj/unit.kg/unit.K
            mu_c = (1-qual)*water.Liquid.mu + qual*water.Vapor.mu
            k_c = (1-qual)*water.Liquid.k + qual*water.Vapor.k
            prtl_c = (1-qual)*water.Liquid.Prandt + qual*water.Vapor.Prandt
        elif water.phase == 'Liquid':
            cp_c = water.cp * unit.kj/unit.kg/unit.K
            mu_c = water.mu
            k_c = water.k
            prtl_c = water.Prandt
        else:
            assert False, 'Vapor not allowed.'

        # Heat transfer coefficient
        diameter = (4*self.core_flow_area/math.pi)**.5
        rey_c = 4*coolant_mass_flowrate / mu_c / math.pi / diameter
        #print('Reactor coolant Rey =',rey_c)
        #print('Reactor coolant mu  =',mu_c)
        #print('Reactor coolant T   =',temp_c-273.15)

        temp_c_w = temp_f
        temp_c_w_F = unit.convert_temperature(temp_c_w, 'K', 'F')
        temp_c_sat_F = unit.convert_temperature(temp_c_sat, 'K', 'F')

        #print('F ',temp_c_w_F, temp_c_sat_F)

        if temp_c_w_F > temp_c_sat_F and \
           coolant_mass_flowrate >= 142.5*unit.kg/unit.second: # nucleate boiling
        # Jens and Lottes correlation for subcooled/saturated nucleate boiling
        # 500 <=  P <= 2000 psi
        # mdot >= 142.5 kg/s

            # Sanity check
            press_c_psia = water.P*unit.mega*unit.pascal/unit.psi
            assert 500 <= press_c_psia <= 2000, 'press_s [psi] = %r'%press_c_psia

            q2prime = ( (temp_c_w_F - temp_c_sat_F) * math.exp(press_c_psia/900) / 60 )**4 * 1e6
            q2prime *= unit.Btu/unit.hour/unit.ft**2
            h_c = q2prime/(temp_c_w - temp_c_sat)
            nusselt_c = h_c * diameter / k_c

            #print('Jens mass flowrate = ',coolant_mass_flowrate)

        else: # single phase heat transfer

            nusselt_c = self.__mean_nusselt_single_phase(rey_c, prtl_c)
            h_c = nusselt_c * k_c / diameter

            #print('Single mass flowrate = ',coolant_mass_flowrate)

        ###################################################
        # Overall heat transfer
        ###################################################

        #therm_cond_wall = self.cladding_k
        # Can we consider the cladding?
        #one_over_U = 1.0/h_c * area_outer/area_inner + \
        #     (radius_outer-radius_inner)/therm_cond_wall * area_outer/area_mean

        one_over_U = 1.0/h_c

        # Total area of heat tranfer
        area = self.core_heat_transfer_area

        UA = area * 1/one_over_U

        # Heating power (W)
        qdot = - UA * (temp_f - temp_c)

        return (qdot, nusselt_c)

    def __mean_nusselt_single_phase(self, rey, prtl):
        """Mean Nusselt number for turbulent one-phase flow on "smooth" pipes.
           Dittus and Boelter.
           Re > 10000
           L/D > 60

           Parameters
           ----------

           rey: float
               Reynolds number based on diameter; bulk conditions
           prtl: float
               Prandtl number; bulk conditions
        """

        if rey <= 2e3:
            nusselt = 4.0
        else:
           coeff = 0.023
           nusselt = coeff * rey**0.8 * prtl**0.4

        return nusselt
