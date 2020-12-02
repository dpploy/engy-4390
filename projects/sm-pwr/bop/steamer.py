#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module.
   Helical coil system steam generator for the NuScale BOP.
   Once-through heat exchanger to operate under fully developed nucleate boiling
   heat transfer. Subcooled? or Saturated? regime??

   + 1012 tubes (Iconel 690)
   + 16 mm OD
   + 0.9 mm tube wall
   + 22.3 m long
   + Primary inflow temperature: 543 F (283.9 C) (w/ a 100 F rise)
   + Primary outflow temperature: 497 F (258.3 C)
   + Primary mass flowrate:
         - Avg: 4.66e+6 (avg) lb/h
         - Max: 5.24e6
         - Min: 4.27e6
   + Secondary inflow temperature: 149 C
   + Secondary outflow temperature: > 241 C
   + Heat transfer area:  17928 ft^2  (1665.57 m^2)
   + Heat flux at operating condition: 8.4077 Btu/ft^2-s (95482.27 W/m^2-s)
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
        self.discard_tau_recording_before = 2*unit.minute
        self.heat_transfer_area = 1665.57*unit.meter**2

        self.helicoil_outer_radius = 16/2*unit.milli*unit.meter
        self.helicoil_tube_wall = 0.9*unit.milli*unit.meter
        self.helicoil_inner_radius = self.helicoil_outer_radius - self.helicoil_tube_wall
        self.helicoil_length = 22.3*unit.meter
        self.n_helicoil_tubes = 3*1012

        self.wall_temp_delta_primary = 1.5*unit.K
        self.wall_temp_delta_secondary = 1.5*unit.K

        self.iconel690_k = 12.1*unit.watt/unit.meter/unit.kelvin

        self.primary_volume = 2.5*unit.meter**3

        self.secondary_volume = math.pi * self.helicoil_inner_radius**2 * \
                                self.helicoil_length * self.n_helicoil_tubes

       # Initialization
        self.primary_inflow_temp = (283.9+273.15)*unit.kelvin
        self.primary_pressure = 127.6*unit.bar
        self.primary_mass_flowrate = 2.5*600*unit.kg/unit.second

        self.primary_outflow_temp = self.primary_inflow_temp #- 2*unit.K

        self.secondary_inflow_temp = (149+273.15)*unit.kelvin

        self.secondary_pressure = 34*unit.bar
        self.secondary_mass_flowrate = 5.4*67*unit.kg/unit.second

        self.secondary_outflow_temp = self.secondary_inflow_temp #- 2*unit.K

        self.secondary_outflow_quality = 0

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
                            value=self.secondary_mass_flowrate,
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
                         formal_name='P_2', unit='Pa',
                         value=self.secondary_pressure,
                         latex_name=r'$P_2$',
                         info='Steamer Secondary Outflow Pressure')

        quantities.append(press)

        quality = Quantity(name='quality',
                         formal_name='X', unit=' ',
                         value=self.secondary_outflow_quality,
                         latex_name=r'$\chi$',
                         info='Steamer Outlet Quality')

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

        self.secondary_outflow_phase.add_row(time, secondary_outflow)
        self.secondary_outflow_phase.set_value('temp', temp_s, time)
        self.secondary_outflow_phase.set_value('flowrate', self.secondary_mass_flowrate, time)
        self.secondary_outflow_phase.set_value('pressure', self.secondary_pressure, time)

        self.state_phase.add_row(time, steamer)

        # Primary properties
        water_p = WaterProps(T=temp_p, P=self.primary_pressure/unit.mega/unit.pascal)
        if water_p.phase == 'Two phases':
            qual = water_p.x
            assert qual <= 0.4 # limit to low quality
            rho_p = (1-qual)*water_p.Liquid.rho + qual*water_p.Vapor.rho
            cp_p = (1-qual)*water_p.Liquid.cp + qual*water_p.Vapor.cp
            cp_p *= unit.kj/unit.kg/unit.K
            mu_p = (1-qual)*water_p.Liquid.mu + qual*water_p.Vapor.mu
            k_p = (1-qual)*water_p.Liquid.k + qual*water_p.Vapor.k
            prtl_p = (1-qual)*water_p.Liquid.Prandt + qual*water_p.Vapor.Prandt
        elif water_p.phase == 'Liquid':
            rho_p = water_p.Liquid.rho
            cp_p = water_p.Liquid.cp*unit.kj/unit.kg/unit.K
            k_p = water_p.Liquid.k
            mu_p = water_p.Liquid.mu
            prtl_p = water_p.Liquid.Prandt
        else:
            assert False,'Vapor not allowed.'

        # Primary residence time
        q_vol_p = self.primary_mass_flowrate/rho_p
        tau_p = self.primary_volume/q_vol_p \
                if q_vol_p > 0 and time > self.discard_tau_recording_before else 0

        self.state_phase.set_value('tau_p', tau_p, time)

        # Secondary properties
        water_s = WaterProps(T=temp_s, P=self.secondary_pressure/unit.mega/unit.pascal)
        if water_s.phase == 'Two phases':
            qual = water.x
            assert qual <= 0.4 # limit to low quality
            rho_s = (1-qual)*water_s.Liquid.rho + qual*water_s.Vapor.rho
            cp_s = (1-qual)*water_s.Liquid.cp + qual*water_s.Vapor.cp
            cp_s *= unit.kj/unit.kg/unit.K
            mu_s = (1-qual)*water_s.Liquid.mu + qual*water_s.Vapor.mu
            k_s = (1-qual)*water_s.Liquid.k + qual*water_s.Vapor.k
            prtl_s = (1-qual)*water_s.Liquid.Prandt + qual*water_s.Vapor.Prandt
        elif water_s.phase == 'Liquid':
            rho_s = water_s.Liquid.rho
            cp_s = water_s.Liquid.cp*unit.kj/unit.kg/unit.K
            k_s = water_s.Liquid.k
            mu_s = water_s.Liquid.mu
            prtl_s = water_s.Liquid.Prandt
        else:
            rho_s = water_s.Vapor.rho
            cp_s = water_s.Vapor.cp*unit.kj/unit.kg/unit.K
            k_s = water_s.Vapor.k
            mu_s = water_s.Vapor.mu
            prtl_s = water_s.Vapor.Prandt

        # Secondary residence time
        q_vol_s = self.secondary_mass_flowrate/rho_s
        tau_s = self.secondary_volume/q_vol_s \
                if q_vol_s > 0 and time > self.discard_tau_recording_before else 0

        self.state_phase.set_value('tau_s', tau_s, time)

        # Heat flux and Nusselt number
        (heat_sink_pwr, nusselt_p, nusselt_s) = self.__heat_sink_pwr(water_p, water_s)

        q2prime = -heat_sink_pwr/self.heat_transfer_area
        if q2prime < 0.0: # case when primary is heated by the secondary
            q2prime = 0.0
        self.state_phase.set_value('heatflux', q2prime, time)

        self.state_phase.set_value('nusselt_p', nusselt_p, time)

        self.state_phase.set_value('nusselt_s', nusselt_s, time)

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
        #print('primary   outflow temp [K] =', temp_p)

        temp_s = u_vec[1] # get temperature of secondary outflow
        #print('secondary outflow temp [K] =', temp_s)

        # initialize f_vec to zero
        f_tmp = np.zeros(2, dtype=np.float64) # vector for f_vec return

        #-----------------------
        # primary energy balance
        #-----------------------

        temp_p_in = self.primary_inflow_temp
        #print('primary inflow temp [K] =', temp_p_in)
        #print('secondary inflow temp [K] =', self.secondary_inflow_temp)

        press_p = self.primary_pressure

        water_p = WaterProps(T=temp_p, P=press_p/unit.mega/unit.pascal)

        assert water_p.phase != 'Two phases'
        assert water_p.phase != 'Vapour'

        rho_p = water_p.rho
        cp_p = water_p.Liquid.cp * unit.kj/unit.kg/unit.K

        vol_p = self.primary_volume
        q_p = self.primary_mass_flowrate/rho_p

        if q_p > 0:
            tau_p = vol_p/q_p
        else:
            tau_p = 1*unit.hour

        #-----------------------
        # secondary energy balance
        #-----------------------

        temp_s_in = self.secondary_inflow_temp
        #print('secondary inflow T [K] =', temp_s_in)

        press_s = self.secondary_pressure
        #print('secondary inflow P [bar] =', press_s/unit.bar)

        #print('secondary P [bar] =', press_s/unit.bar)
        #print('Primary   T [K]   =', temp_p)
        #print('Secondary T [K]   =', temp_s)
        #print('quality = ', self.secondary_outflow_quality)

        water_s = WaterProps(T=temp_s, P=press_s/unit.mega/unit.pascal)

        if water_s.phase == 'Liquid':
            rho_s = water_s.Liquid.rho
            cp_s = water_s.Liquid.cp
        elif water_s.phase == 'Two phases':
            cp_s = (1-qual)*water.Liquid.cp + qual*water.Vapor.cp
            rho_s = (1-qual)*water.Liquid.rho + qual*water.Vapor.rho
        else:
            rho_s = water_s.Vapor.rho
            cp_s = water_s.Vapor.cp

        cp_s *= unit.kj/unit.kg/unit.K

        vol_s = self.secondary_volume
        q_s = self.secondary_mass_flowrate/rho_s

        if q_s > 0:
            tau_s = vol_s/q_s
        else:
            tau_s = 1*unit.hour

        #-----------------------
        # calculations
        #-----------------------
        (heat_sink_pwr, _, _) = self.__heat_sink_pwr(water_p, water_s)
        heat_sink_pwr_dens = heat_sink_pwr/vol_p

        #print('heat_sink =', heat_sink_pwr_dens)

        #assert heat_sink_pwr < 0, 'heatsink =  %r'%(heat_sink_pwr)
        #assert temp_p-temp_p_in < 0

        f_tmp[0] = - 1/tau_p * (temp_p - temp_p_in) +  heat_sink_pwr_dens/(rho_p*cp_p)

        heat_source_pwr = - heat_sink_pwr
        heat_source_pwr_dens = heat_source_pwr/vol_s

        #assert temp_s-temp_s_in < 0
        f_tmp[1] = - 1/tau_s * (temp_s - temp_s_in) + heat_source_pwr_dens/(rho_s*cp_s)

        #print('Primary:',f_tmp[0], 'Convective = ', - 1/tau_p * (temp_p - temp_p_in), 'heat term = ', heat_sink, 'mass*cp = ' ,(rho_p*cp_p*vol_p) )
        #print('Secondary:',f_tmp[1], 'Convective = ', - 1/tau_s * (temp_s - temp_s_in), 'heat term = ', heat_source ,'mass*cp = ', (rho_s*cp_s*vol_s))


        # Boiling test

        sat_liq = WaterProps(P=press_s/unit.mega/unit.pascal, x=0.0)
        sat_vap = WaterProps(P=press_s/unit.mega/unit.pascal, x=1.0)
        sensible_water = WaterProps(T= (temp_s_in +sat_liq.T)/2 , P=press_s/unit.mega/unit.pascal)
        cp_l_o = sensible_water.cp
        q_total = (temp_s - temp_s_in)*cp_s/1000
        #q_total = heat_source/1000/self.secondary_mass_flowrate
        q_heat = (sat_liq.T- temp_s_in)*cp_l_o
        h_v = sat_vap.h
        h_l = sat_liq.h
        h_vap = h_v-h_l
        q_vap = h_vap
        if q_total < q_heat: # subcooled
            self.secondary_outflow_quality = 0
        elif q_total > (q_vap+q_heat): #superheated
            self.secondary_outflow_quality = 1
        else: #mixed
            self.secondary_outflow_quality = 1- (h_vap - (q_total - q_heat))/h_vap

        return f_tmp

    def __heat_sink_pwr(self, water_p, water_s):
        """Cooling rate of the primary side.

           Assumptions
           -----------

           + primary side: overall single phase heat tranfer. Locally there may be
             either partial nucleate boiling or fully developed nucleate boiling
             but the model will not capture this.

           + secondary side: overall ranging from one phase heat tranfer to fully
             developed nucleate boiling .
        """
        # Primary props
        temp_p = water_p.T

        water_p_sat = WaterProps(P=water_p.P, x=0.0)
        temp_p_sat = water_p_sat.T

        # Overall condition on the primary; locally there may be nucleate boiling
        assert temp_p < temp_p_sat

        cp_p = water_p.cp
        mu_p = water_p.mu
        k_p = water_p.k
        prtl_p = water_p.Prandt
        cp_p *= unit.kj/unit.kg/unit.K

        # Secondary props
        temp_s = water_s.T

        water_s_sat = WaterProps(P=water_s.P, x=0.0)
        temp_s_sat = water_s_sat.T

        # Overall condition on the secondary
        #assert temp_s <= temp_s_sat + 80*unit.kelvin # CHF calc neeeded

        if water_s.phase == 'Two phases':
            qual = water_s.x
            cp_s = (1-qual)*water_s.Liquid.cp + qual*water_s.Vapor.cp
            mu_s = (1-qual)*water_s.Liquid.mu + qual*water_s.Vapor.mu
            k_s = (1-qual)*water_s.Liquid.k + qual*water_s.Vapor.k
            prtl_s = (1-qual)*water_s.Liquid.Prandt + qual*water_s.Vapor.Prandt
        else:
            cp_s = water_s.cp
            mu_s = water_s.mu
            k_s = water_s.k
            prtl_s = water_s.Prandt

        cp_s *= unit.kj/unit.kg/unit.K

        ###########################################
        # Heat transfer coefficient on primary side
        ###########################################
        radius_outer = self.helicoil_outer_radius

        rey_p = self.primary_mass_flowrate * 2*radius_outer / mu_p
        #print('primary Reynolds = ',rey_p)

        sl = 1.0      # tube bundle pitch parallel to flow
        st = 1.5 * sl # tube bundle pitch transverse to flow

        temp_p_w = self.primary_inflow_temp - self.wall_temp_delta_primary # wall temperature
        #print('Prinary Wall Temp [K] =',temp_p_w)

        water_p_w = WaterProps(T=temp_p_w, P=water_p.P) # primary at wall T, P

        assert water_p_w.phase != 'Two phases' # sanity check

        prtl_w = water_p_w.Prandt

        nusselt_p = self.__mean_nusselt_single_phase(rey_p, prtl_p, prtl_w, st/sl)

        h_p = nusselt_p * k_p / (2*radius_outer)

        #############################################
        # Heat transfer coefficient on secondary side
        #############################################
        radius_inner = self.helicoil_inner_radius

        #temp_s_w = temp_p_w - self.wall_temp_delta_secondary
        temp_s_w = temp_p_w

        water_s_sat = WaterProps(P=water_s.P, x=0.0)
        temp_s_sat = water_s_sat.T

        #print('secondary pressure [bar] = ',water_s.P*unit.mega*unit.pascal/unit.bar)
        #print('primary   pressure [bar] = ',water_p.P*unit.mega*unit.pascal/unit.bar)

        #assert temp_s_w >= temp_s, 'temp_s = %r, temp_s_w = %r'%(temp_s,temp_s_w)
        #print('saturated temp = ', temp_s_sat)

        if (temp_s_w - temp_s_sat) > 0.0: # nucleate boiling
        # Jens and Lottes correlation for subcooled/saturated nucleate boiling
        # 3.5 <= P <= 14 MPa
            q2prime = ((temp_s_w - temp_s_sat)*math.exp(self.secondary_pressure/unit.mega/unit.pascal/6.2)/0.79)**4
            h_s = q2prime/(temp_s_w - temp_s_sat)
            nusselt_s = h_s * 2*radius_inner / k_s
        else: # single phase transfer
            rey_s = self.secondary_mass_flowrate * 2*radius_inner / mu_s
            water_s_w = WaterProps(T=temp_s_w, P=water_s.P) # secondary at wall T, P
            assert water_s_w.phase != 'Two phases' # sanity check
            prtl_w = water_s_w.Prandt
            nusselt_s = self.__mean_nusselt_single_phase(rey_s, prtl_s, prtl_w, st/sl)
            h_s = nusselt_s * k_s / (2*radius_inner)

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
        if h_s > 0:
            one_over_U = 1.0/h_p * area_outer/area_inner + \
            (radius_outer-radius_inner)/therm_cond_wall * area_outer/area_mean + \
            fouling + \
            1/h_s
        else:
            one_over_U = 1.0/h_p * area_outer/area_inner + \
            (radius_outer-radius_inner)/therm_cond_wall * area_outer/area_mean + \
            fouling

        # Total area of heat tranfer
        #area = 2*math.pi*radius_mean* self.n_helicoil_tubes * self.helicoil_length
        area = self.heat_transfer_area

        temp_p_avg = (self.primary_inflow_temp + temp_p)/2
        temp_s_avg = (self.secondary_inflow_temp + temp_s)/2
        qdot = - area * 1/one_over_U * (temp_p_avg-temp_s_avg)
        qdot = - 95482.27 * 1665.57
        #qdot = - area * 1/one_over_U * (self.primary_inflow_temp-self.secondary_inflow_temp)
        #print('here',qdot/area)

        return (qdot, nusselt_s, nusselt_p)

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
