#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module.
   Helical coil system steam generator for the NuScale BOP.
   Once-through heat exchanger to operate under fully developed nucleate boiling
   heat transfer. Subcooled? or Saturated? regime??

   + 1380 tubes (Iconel 690)
   + 16 mm OD
   + 0.9 mm tube wall
   + 22.3 m long
   + Primary inflow temperature: 543 F (283.9 C) (w/ a 100 F rise)
   + Primary outflow temperature: 497 F (258.3 C)
   + Primary mass flowrate:
         - Avg: 4.66e+6 (avg) lb/h (587.15 kg/s)
         - Max: 5.24e6 (600.22 kg/s)
         - Min: 4.27e6 (538.01 kg/s)
   + Secondary inflow temperature: 149 C
   + Secondary outflow temperature: 584.4 F (306.9 C) (saturated = 241.68 C)
   + Heat transfer area:  17928 ft^2  (1665.57 m^2)
   + Heat flux at operating condition: 8.4077 Btu/ft^2-s (95482.27 W/m^2)
   + Full stem flow: 532100 lb/h (67.04 kg/s)
   + Operating secondary temperature: 575 F (301.67 C)
   + Operating secondary pressure: 500 psi (34.47 bar)

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

class Steamer(Module):
    """Steam generator.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: reactor, turbine.
    See instance attribute `port_names_expected`.

    """

    def __init__(self, primary_inflow_temp=20+273.15, secondary_inflow_temp=20+273.15):
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

        # Initialization
        self.primary_inflow_temp = primary_inflow_temp

        self.primary_pressure = 127.6*unit.bar

        self.primary_mass_flowrate = 4.66e6*unit.lb/unit.hour

        self.primary_outflow_temp = self.primary_inflow_temp #- 2*unit.K

        self.secondary_inflow_temp = secondary_inflow_temp

        self.secondary_pressure = 34*unit.bar
        self.secondary_mass_flowrate = 67*unit.kg/unit.second

        self.secondary_outflow_temp = self.secondary_inflow_temp #- 2*unit.K

        self.secondary_outflow_quality = 0 # running value of quality

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

        # Primary outflow phase history
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

    def __step(self, time=0.0):
        """ODE IVP problem.
        """

        # Get state values
        u_0 = self.__get_state_vector(time)

        t_interval_sec = np.linspace(time, time+self.time_step, num=2)

        max_n_steps_per_time_step = 1500 # max number of nonlinear algebraic solver
                                         # iterations per time step

        (u_vec_hist, info_dict) = odeint(self.__f_vec, u_0, t_interval_sec,
                                         #rtol=1e-4, atol=1e-8,
                                         mxstep=max_n_steps_per_time_step,
                                         full_output=True, tfirst=False)

        assert info_dict['message'] == 'Integration successful.', info_dict['message']

        u_vec = u_vec_hist[1, :]  # solution vector at final time step

        temp_p = u_vec[0] # primary outflow temp
        temp_s = u_vec[1] # secondary outflow temp

        # Update phases
        primary_outflow = self.primary_outflow_phase.get_row(time)
        secondary_outflow = self.secondary_outflow_phase.get_row(time)
        steamer = self.state_phase.get_row(time)

        time += self.time_step

        self.primary_outflow_phase.add_row(time, primary_outflow)
        self.primary_outflow_phase.set_value('temp', temp_p, time)
        self.primary_outflow_phase.set_value('flowrate', self.primary_mass_flowrate, time)

        self.secondary_outflow_phase.add_row(time, secondary_outflow)
        self.secondary_outflow_phase.set_value('temp', temp_s, time)
        self.secondary_outflow_phase.set_value('flowrate', self.secondary_mass_flowrate, time)
        self.secondary_outflow_phase.set_value('pressure', self.secondary_pressure, time)
        self.secondary_outflow_phase.set_value('quality', self.secondary_outflow_quality, time)

        self.state_phase.add_row(time, steamer)

        # Primary residence time
        self.state_phase.set_value('tau_p', self.tau_p, time)

        # Secondary residence time
        self.state_phase.set_value('tau_s', self.tau_s, time)

        # Heat flux and Nusselt number
        heatflux = -self.heat_sink_pwr/self.heat_transfer_area
        self.state_phase.set_value('heatflux', heatflux, time)

        self.state_phase.set_value('nusselt_p', self.nusselt_p, time)

        self.state_phase.set_value('nusselt_s', self.nusselt_s, time)

        return time

    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.
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

        # Heat transfered
        self.heat_sink_pwr = self.__heat_sink_pwr(temp_p, temp_s)

        #-----------------------
        # primary energy balance
        #-----------------------

        temp_p_in = self.primary_inflow_temp

        vol_p = self.primary_volume
        q_p = self.primary_mass_flowrate/self.rho_p

        if q_p > 0:
            tau_p = vol_p/q_p
        else:
            tau_p = 1*unit.hour

        self.tau_p = tau_p

        #-------------------------
        # Secondary energy balance
        #-------------------------

        temp_s_in = self.secondary_inflow_temp

        vol_s = self.secondary_volume
        q_s = self.secondary_mass_flowrate/self.rho_s

        if q_s > 0:
            if self.secondary_outflow_quality == 0:
                tau_s = vol_s/q_s  # liquid residence time
            elif self.secondary_outflow_quality == 1:
                #tau_s = 0.85*vol_s/q_s # faster moving gas with higher flow losses
                # this controls how high temp_s will jumpt to
                tau_s = 0.35*vol_s/q_s # vapor flow acceleration
            else: # vapor/liquid mix
                # This factor controls the onset of the transition jump in temp_s
                tau_s = 0.9*vol_s/q_s # vapor/liquid flow acceleration
        else:
            tau_s = 1*unit.hour

        self.tau_s = tau_s

        #----------------------------
        # Compute ODE right-hand side
        #----------------------------

        heat_sink_pwr_dens = self.heat_sink_pwr/vol_p

        f_tmp[0] = - 1/tau_p * (temp_p - temp_p_in) + \
                   heat_sink_pwr_dens/(self.rho_p*self.cp_p)

        heat_source_pwr = - self.heat_sink_pwr

        heat_source_pwr -= self.__heat_vaporization_pwr()

        heat_source_pwr_dens = heat_source_pwr/vol_s

        f_tmp[1] = - 1/tau_s * (temp_s - temp_s_in) + \
                   heat_source_pwr_dens/(self.rho_s*self.cp_s)

        #-----------------------
        # Wrap up
        #-----------------------

        self.__update_quality(temp_p, temp_s)

        return f_tmp

    def __heat_sink_pwr(self, temp_p, temp_s):
        """Cooling rate of the primary side.

           Assumptions
           -----------

           + primary side: overall single phase heat transfer. Locally there may be
             either partial nucleate boiling or fully developed nucleate boiling
             but the model will not capture this.

           + secondary side: overall ranging from one phase heat transfer to fully
             developed nucleate boiling.
        """

        ###########################################
        # Heat transfer coefficient on primary side
        ###########################################

        h_p = self.__heat_transfer_coeff_primary(temp_p)

        #############################################
        # Heat transfer coefficient on secondary side
        #############################################

        h_s = self.__heat_transfer_coeff_secondary(temp_p, temp_s)

        ###################################################
        # Overall heat transfer
        ###################################################
        # This is based on the secondary side

        # Cross section geometry
        area_outer = math.pi * self.helicoil_outer_radius**2
        area_inner = math.pi * self.helicoil_inner_radius**2
        radius_mean = (self.helicoil_outer_radius+self.helicoil_inner_radius)/2
        area_mean = math.pi * radius_mean**2
        radius_inner = self.helicoil_inner_radius
        radius_outer = self.helicoil_outer_radius
        therm_cond_wall = self.iconel690_k

        fouling = 0.0003 * unit.F*unit.ft**2*unit.hour/unit.Btu

        # Secondary side based heat transfer resistance
        if h_s > 1:
            one_over_U = 1.0/h_p * area_outer/area_inner + \
            (radius_outer-radius_inner)/therm_cond_wall * area_outer/area_mean + \
            fouling + \
            1/h_s
        else:
            one_over_U = 1.0/h_p * area_outer/area_inner + \
            (radius_outer-radius_inner)/therm_cond_wall * area_outer/area_mean + \
            fouling

        # Total area of heat transfer
        #area = 2*math.pi*radius_mean* self.n_helicoil_tubes * self.helicoil_length
        area = self.heat_transfer_area

        temp_p_avg = (self.primary_inflow_temp + temp_p)/2
        temp_s_avg = (self.secondary_inflow_temp + temp_s)/2

        qdot = - area * 1/one_over_U * (temp_p_avg - temp_s_avg)

        qdot /= 2.61 # calibration factor for correlation to work with NuScale specs

        #print(qdot, -95482.27 * 1665.57)
        #qdot = - 95482.27 * 1665.57

        return qdot

    def __heat_transfer_coeff_primary(self, temp_p):
        """Heat transfer coefficient in the primary flow. SI units: W/m2-K.
        """

        # Primary props
        press_p_MPa = self.primary_pressure/unit.mega/unit.pascal
        water_p = WaterProps(T=temp_p, P=press_p_MPa)

        assert water_p.phase == 'Liquid' # sanity check

        self.rho_p = water_p.rho
        self.cp_p = water_p.Liquid.cp * unit.kj/unit.kg/unit.K
        self.mu_p = water_p.Liquid.mu
        self.k_p = water_p.Liquid.k
        self.prtl_p = water_p.Liquid.Prandt

        radius_outer = self.helicoil_outer_radius

        self.rey_p = self.primary_mass_flowrate * 2*radius_outer / self.mu_p

        temp_p_w = temp_p - self.wall_temp_delta_primary # wall temperature

        water_p_w = WaterProps(T=temp_p_w, P=press_p_MPa) # @ wall

        assert water_p_w.phase == 'Liquid' # sanity check

        prtl_w = water_p_w.Liquid.Prandt

        st_over_sl = self.tube_bundle_pitch_ratio

        self.nusselt_p = self.__mean_nusselt_single_phase(self.rey_p, self.prtl_p, prtl_w, st_over_sl)

        h_p = self.nusselt_p * self.k_p / (2*radius_outer)

        return h_p

    def __mean_nusselt_single_phase(self, rey, prtl, prtl_w, st_sl):
        """Mean Nusselt number for turbulent one-phase flow.
           Staggered tube bundle.
           A. Zukauskas 1987,
           Convective Heat Transfer in Cross Flows
           Handbook of Single-Phase Convective Heat Transfer, Chap 6.
           S. Kakac, R. Shah, and W. Aung Eds.
           J. Wiley & Sons, New York 1987

           Parameters
           ----------

           rey: float
               Reynolds number based on diameter
           prtl: float
               Prandtl number
           prtl_w: float
               Prandtl number based on wall temperature.
           st_sl: float
               Ratio of tube pitch transversal to flow to the tube pitch along flow.
        """

        if rey < 1:
            nusselt_p = 4.0
        elif 1 <= rey <= 5e2:
            nusselt_p = 1.04 * rey**0.4 * prtl**0.36 * \
                        (prtl/prtl_w)**0.25
        elif 5e2 < rey <= 1e3:
            nusselt_p = 0.71 * rey**0.5 * prtl**0.36 * \
                        (prtl/prtl_w)**0.25
        elif 1e3 < rey <= 2e5:
            nusselt_p = 0.35 * (st_sl)**0.2 * rey**0.6 * prtl**0.36 * \
                        (prtl/prtl_w)**0.25
        elif 2e5 < rey <= 2e6:
            nusselt_p = 0.031 * (st_sl)**0.2 * rey**0.8 * prtl**0.36 * \
                        (prtl/prtl_w)**0.25
        else:
            assert False

        return nusselt_p

    def __heat_transfer_coeff_secondary(self, temp_p, temp_s):
        """Heat transfer coefficient in the secondary flow. SI units: W/m2-K.
        """

        press_s_MPa = self.secondary_pressure/unit.mega/unit.pascal

        if self.secondary_outflow_quality == 0:
            water_s = WaterProps(T=temp_s, P=press_s_MPa)
            self.rho_s = water_s.rho
            cp_s = water_s.cp
            self.mu_s = water_s.mu
            self.k_s = water_s.k
            self.prtl_s = water_s.Prandt
        elif self.secondary_outflow_quality == 1:
            water_s = WaterProps(T=temp_s, P=press_s_MPa)
            self.rho_s = water_s.rho
            cp_s = water_s.cp
            self.mu_s = water_s.mu
            self.k_s = water_s.k
            self.prtl_s = water_s.Prandt
        else:
            qual = self.secondary_outflow_quality
            w_sat = WaterProps(P=press_s_MPa, x=qual)
            #mass_frac_v = w_sat.Vapor.rho / (w_sat.Liquid.rho + w_sat.Vapor.rho)
            #self.rho_s = (1-mass_frac_v)*w_sat.Liquid.rho + mass_frac_v*w_sat.Vapor.rho
            #cp_s = (1-mass_frac_v)*w_sat.Liquid.cp + mass_frac_v*w_sat.Vapor.cp
            #self.mu_s = (1-mass_frac_v)*w_sat.Liquid.mu + mass_frac_v*w_sat.Vapor.mu
            #self.k_s = (1-mass_frac_v)*w_sat.Liquid.k + mass_frac_v*w_sat.Vapor.k

            self.rho_s = w_sat.rho
            cp_s = (1-qual)*w_sat.Liquid.cp + qual*w_sat.Vapor.cp
            self.mu_s = (1-qual)*w_sat.Liquid.mu + qual*w_sat.Vapor.mu
            self.k_s = (1-qual)*w_sat.Liquid.k + qual*w_sat.Vapor.k

            self.prtl_s = cp_s*unit.kj/unit.kg/unit.K*self.mu_s/self.k_s

        cp_s *= unit.kj/unit.kg/unit.K
        self.cp_s = cp_s

        radius_inner = self.helicoil_inner_radius

        self.rey_s = self.secondary_mass_flowrate * 2*radius_inner / self.mu_s

        water_s_sat = WaterProps(P=press_s_MPa, x=0.5)
        temp_s_sat = water_s_sat.T

        # Overall condition on the secondary
        #assert temp_s <= temp_s_sat + 80*unit.kelvin # CHF calc neeeded

        temp_p_w = temp_p - self.wall_temp_delta_primary # wall temperature
        temp_s_w = temp_p_w - self.wall_temp_delta_secondary

        #assert temp_s_w >= temp_s, 'temp_s = %r, temp_s_w = %r'%(temp_s,temp_s_w)
        #print('saturated temp = ', temp_s_sat)

        if (temp_s_w - temp_s_sat) > 5: # nucleate boiling (allow for development)
        # Jens and Lottes correlation for subcooled/saturated nucleate boiling
        # 3.5 <= P <= 14 MPa
            q2prime = ((temp_s_w - temp_s_sat)*math.exp(self.secondary_pressure/unit.mega/unit.pascal/6.2)/0.79)**4
            h_s = q2prime/(temp_s_w - temp_s_sat)
            self.nusselt_s = h_s * 2*radius_inner / self.k_s
            #print('q2prime=',q2prime/unit.kilo,'DT=',temp_s_w-temp_s_sat)
            #print('h_s JL=',h_s)
        else: # single phase transfer
            #print(temp_s_w-273.15, press_s_MPa)
            water_s_w = WaterProps(T=temp_s_w, P=press_s_MPa) # @ wall
            assert water_s_w.phase != 'Two phases' # sanity check
            if water_s_w.phase == 'Liquid':
                prtl_w = water_s_w.Liquid.Prandt
            elif water_s_w.phase == 'Vapour':
                prtl_w = water_s_w.Vapor.Prandt
            else:
                assert False
            st_over_sl = self.tube_bundle_pitch_ratio
            self.nusselt_s = self.__mean_nusselt_single_phase(self.rey_s, self.prtl_s, prtl_w, st_over_sl)
            h_s = self.nusselt_s * self.k_s / (2*radius_inner)

        return h_s

    def __heat_vaporization_pwr(self):
        """Heat of vaporization in the secondary.
        """

        press_p_MPa = self.primary_pressure/unit.mega/unit.pascal

        sat_liq = WaterProps(P=press_p_MPa, x=0.0)
        sat_vap = WaterProps(P=press_p_MPa, x=1.0)

        h_v = sat_vap.h
        h_l = sat_liq.h
        h_vap = h_v - h_l

        return self.secondary_outflow_quality * h_vap*unit.kj * self.secondary_mass_flowrate

    def __update_quality(self, temp_p, temp_s):

        press_s_MPa = self.secondary_pressure/unit.mega/unit.pascal

        # Secondary steam quality
        sat_liq = WaterProps(P=press_s_MPa, x=0.0)
        sat     = WaterProps(P=press_s_MPa, x=0.5)
        sat_vap = WaterProps(P=press_s_MPa, x=1.0)

        temp_s_in = self.secondary_inflow_temp

        assert temp_s_in < sat.T
        sensible_water = WaterProps(T=(temp_s_in + sat_liq.T)/2 , P=press_s_MPa)
        cp_sensible = sensible_water.Liquid.cp*unit.kj/unit.kg/unit.K

        q_sensible = (sat_liq.T - temp_s_in)*cp_sensible

        spcf_heat_transfered = - self.heat_sink_pwr/self.secondary_mass_flowrate
        #print("q''[kW/m2]=",-self.heat_sink_pwr/self.heat_transfer_area/unit.kilo)

        h_v = sat_vap.h * unit.kj/unit.kg
        h_l = sat_liq.h * unit.kj/unit.kg
        h_vap = h_v - h_l

        if spcf_heat_transfered < q_sensible and temp_s < sat.T: # subcooled
            self.secondary_outflow_quality = 0
        elif spcf_heat_transfered > q_sensible + h_vap and temp_s > sat.T: # superheated
            self.secondary_outflow_quality = 1
        elif spcf_heat_transfered > q_sensible + h_vap and temp_s <= sat.T: # slowed superheated
            # Heat transfered is sufficient in equilibrium but in non-equilibrium,
            # there is inertia.
            qual = temp_s/sat.T
            self.secondary_outflow_quality = min(qual,0.999)
            qual = self.secondary_outflow_quality
            w_sat = WaterProps(P=press_s_MPa, x=qual)
            #mass_frac_v = w_sat.Vapor.rho / (w_sat.Liquid.rho + w_sat.Vapor.rho)
            #cp_s = (1-mass_frac_v)*w_sat.Liquid.cp + mass_frac_v*w_sat.Vapor.cp
            cp_s = (1-qual)*w_sat.Liquid.cp + qual*w_sat.Vapor.cp
            cp_s *= unit.kj/unit.kg/unit.K
            #rho_s = (1-mass_frac_v)*w_sat.Liquid.rho + mass_frac_v*w_sat.Vapor.rho
            rho_s = w_sat.rho

        elif q_sensible <= spcf_heat_transfered <= q_sensible + h_vap:
            self.secondary_outflow_quality = (spcf_heat_transfered - q_sensible)/ h_vap
            qual = self.secondary_outflow_quality
            w_sat = WaterProps(P=press_s_MPa, x=qual)
            #mass_frac_v = w_sat.Vapor.rho / (w_sat.Liquid.rho + w_sat.Vapor.rho)
            #cp_s = (1-mass_frac_v)*w_sat.Liquid.cp + mass_frac_v*w_sat.Vapor.cp
            cp_s = (1-qual)*w_sat.Liquid.cp + qual*w_sat.Vapor.cp
            cp_s *= unit.kj/unit.kg
            #rho_s = (1-mass_frac_v)*w_sat.Liquid.rho + mass_frac_v*w_sat.Vapor.rho
            rho_s = w_sat.rho
        else:
            print('ht %r, q_sens %r, h_vap %r, temp_s %r, temp_sat %r'%(spcf_heat_transfered,q_sensible,h_vap,temp_s,sat.T))
            assert False, 'bailing out; unknown case.'
