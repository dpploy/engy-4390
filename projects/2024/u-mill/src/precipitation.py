#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""
Cortix Module
This module is a model of the Precipitation process in the White Mesa Uranium Milling Plant Project


           Feed (Uranyl Tri-Sulfate Stripping product from SolvEx)
             |
             |
             V
   ------------------------
   |                      |
   |    Precipitation     |<------ NH4 (anhydrous) + air
   |                      |
   |----------------------|
   |                      |
   | Thickening / Washing |<------ wash water centrifuge
   |                      |
   ------------------------
             |
             |
             |
             v
         Ammonium
     Diuranate Product


 + Precipitation
   Continuous precipitation with anhydrous ammonia.
   Receives UO2(SO4)3^4- solution from SolvEx; sparge with 1:3 NH4/Air; kept at 30 to 50 C and
   terminal 7 <= pH <= 8:.

       2 UO2(SO4)3^4- + 6 NH3 -> (NH4)2U2O7 + 4 SO4^2-

   Ammonium diuranate (NH4)2U2O7) has a diuranate anion U2O7^2-.

 + Thickening (2 solid bowl centrifuges in series)
   - Wash water volumetric flowrate


Source:
 1993 Uranium Extraction Technology, IAEA Technical Reports Series No. 359
  p. 338 (White Mesa),
  p. 240 (Precipitation: the ammonia system)
 URL: https://www-pub.iaea.org/MTCD/Publications/PDF/trs359_web.pdf

"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity
from cortix import Species

from cortix import Units as unit

class Precipitation(Module):
    """Precipitation.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: Solvex, Evaporator.
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
    #orig def __init__(self, primary_inflow_temp=20+273.15, secondary_inflow_temp=20+273.15):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['uts-feed', 'adu-product']

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)
        self.save = True

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters
        self.precipitation_tank_vol = 45 * unit.meter**3
        self.thickening_tank_vol = 45 * unit.meter**3

        # Initialization

        # Precipitation
        self.uranyl_trisulfate_feed_mass_flowrate = 727 * unit.kg/unit.minute
        self.uranyl_trisulfate_feed_mass_density = 1.1 * unit.kg / unit.liter

        self.ammonia_sparge_feed_mass_flowrate = 5 * unit.kg/unit.minute
        self.ammonia_sparge_feed_mass_density = 5 * unit.kg/unit.minute

        # Thickening/Washing

        '''
        self.primary_inflow_temp = primary_inflow_temp

        self.primary_pressure = 127.6*unit.bar

        self.primary_mass_flowrate = 4.66e6*unit.lb/unit.hour

        self.primary_outflow_temp = self.primary_inflow_temp #- 2*unit.K

        self.secondary_inflow_temp = secondary_inflow_temp

        self.secondary_pressure = 34*unit.bar
        self.secondary_mass_flowrate = 67*unit.kg/unit.second

        self.secondary_outflow_temp = self.secondary_inflow_temp #- 2*unit.K

        self.secondary_outflow_quality = 0 # running value of quality
        '''

        # Derived quantities
        '''
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
        '''

        #***************************************************************************************
        # P R E C I P I T A T I O N
        #***************************************************************************************

        # Precipitation Phase History
        quantities = list()
        species = list()

        precipitation_mass_flowrate = Quantity(name='mass-flowrate',
                                               formal_name='mdot', unit='kg/s',
                                               value=0.0,
                                               latex_name=r'$\dot{m}_{p}$',
                                               info='Precipitation Mass Flowrate')
        quantities.append(precipitation_mass_flowrate)

        self.precipitation_phase = Phase(time_stamp=self.initial_time,
                                         time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # T H I C K N E N I N G
        #***************************************************************************************

        # Thickening/Centrifuge Phase History (thickening/centrifuge)
        quantities = list()
        species = list()

        thickening_mass_flowrate = Quantity(name='mass-flowrate',
                                       formal_name='mdot', unit='kg/s',
                                       value=0.0,
                                       latex_name=r'$\dot{m}_{t}$',
                                       info='Thickening Mass Flowrate')
        quantities.append(thickening_mass_flowrate)

        self.thickening_phase = Phase(time_stamp=self.initial_time,
                                      time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # S T A T E  P H A S E
        #***************************************************************************************
# State Phase History
        quantities = list()
        species = list()

        state_mass_flowrate = Quantity(name='mass-flowrate',
                                        formal_name='mdot', unit='kg/s',
                                        value=0.0,
                                        latex_name=r'$\dot{m}_{s}$',
                                        info='State Mass Flowrate')
        quantities.append(state_mass_flowrate)

        residence_time = Quantity(name='tau_p',
                                   formal_name='tau', unit='s',
                                   value=0.0,
                                   latex_name=r'$\tau_{p}$',
                                   info='Residence Time in Precipitation')
        quantities.append(residence_time)

        heat_flux = Quantity(name='heatflux',
                              formal_name='q', unit='W',
                              value=0.0,
                              latex_name=r'$q$',
                              info='Heat Transfer Rate')
        quantities.append(heat_flux)

        nusselt_number = Quantity(name='nusselt_number',
                                   formal_name='Nu', unit='-',
                                   value=0.0,
                                   latex_name=r'Nu',
                                   info='Nusselt Number')
        quantities.append(nusselt_number)

        self.state_phase = Phase(time_stamp=self.initial_time,
                                 time_unit='s', quantities=quantities, species=species)

    def run(self, *args):

        # Some logic for logging time stamps
        # Leave this here: rebuild logger
        logger_name = args[0][0].name
        self.rebuild_logger(logger_name)

        self.end_time = max(self.end_time, self.initial_time + self.time_step)

        time = self.initial_time

        print_time = self.initial_time
        print_time_step = self.show_time[1]

        print_time_step = max(print_time_step, self.time_step)

        while time <= self.end_time:

            if self.show_time[0] and \
               (print_time <= time < print_time+print_time_step):

                msg = self.name+'::run():time[d]='+ str(round(time/unit.day, 1))
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

        self.end_time = time # correct the final time if needed

    def __call_ports(self, time):

        # Interactions in the uts-feed port
        #----------------------------------
        # One way "from" uts-feed

        # Receive from
        if self.get_port('uts-feed').connected_port:

            self.send(time, 'uts-feed')

            (check_time, uts_feed) = self.recv('uts-feed')
            assert abs(check_time-time) <= 1e-6

            '''
            self.primary_inflow_temp = primary_inflow['temperature']
            self.primary_ressure = primary_inflow['pressure']
            self.primary_mass_flowrate = primary_inflow['mass_flowrate']
            '''

        # Interactions in the adu-port
        #-----------------------------------------
        # One way "to" adu-port

        # Send to
        if self.get_port('adu-product').connected_port:

            msg_time = self.recv('adu-product')

            #temp = self.primary_outflow_phase.get_value('temp', msg_time)

            product = dict()
            '''
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_pressure
            primary_outflow['mass_flowrate'] = self.primary_mass_flowrate
            primary_outflow['quality'] = 0.0
            '''

            self.send((msg_time, product), 'adu-product')

    def __step(self, time=0.0):
        """Stepping Decantation-Filtration in time
        """

        '''
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

        temp_p = u_vec[0] # primary outflow temp
        temp_s = u_vec[1] # secondary outflow temp

        # Update phases
        primary_outflow = self.primary_outflow_phase.get_row(time)
        secondary_inflow = self.secondary_inflow_phase.get_row(time)
        secondary_outflow = self.secondary_outflow_phase.get_row(time)
        steamer = self.state_phase.get_row(time)
        '''
        # Evolve the precipitation state
        mass_flowrate_initial = self.precipitation_phase.get_value('mass-flowrate', time)

        uts_feed_mass_flowrate = self.uranyl_trisulfate_feed_mass_flowrate
        sparge_feed_mass_flowrate = self.ammonia_sparge_feed_mass_flowrate

        # Ideal solution
        mass_flowrate_inflow = uts_feed_mass_flowrate + sparge_feed_mass_flowrate

        rho_uts_feed = self.uranyl_trisulfate_feed_mass_density
        rho_nh4_feed = self.ammonia_sparge_feed_mass_density

        # Ideal solution
        rho_precipitation = rho_uts_feed + rho_nh4_feed

        vol_flowrate_initial = mass_flowrate_initial/rho_precipitation

        if vol_flowrate_initial == 0:
            vol_flowrate_initial = mass_flowrate_inflow/rho_precipitation
            tau = self.precipitation_tank_vol/vol_flowrate_initial
        else:
            tau = self.precipitation_tank_vol/vol_flowrate_initial

        # Mass balance
        mass_flowrate_precipitation = mass_flowrate_inflow + \
                                 math.exp(-time/tau) * (mass_flowrate_initial - mass_flowrate_inflow)

        tmp_precipitation = self.precipitation_phase.get_row(time)

        # Advance time and store new state variables
        time += self.time_step

        self.precipitation_phase.add_row(time, tmp_precipitation)
        self.precipitation_phase.set_value('mass-flowrate', mass_flowrate_precipitation, time)

        '''
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
        '''

        return time
