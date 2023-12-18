#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""
Cortix Module
This module is a model of the Solvent Extraction process in the White Mesa Uranium Milling Plant

   Stripping Feed
(from Precipitation)
            |             |
            |             |  Extraction Feed (from Decantation-Filtration)
            |             |
            V             v
           |----------------|
           |    Solvent     |--------> Raffinate Stream (to Decantation-Filtration)
           |   Extraction   |
           |                |
           |   Scrubbing    |
           |                |<-------- Organic Feed (internal)
           |   Stripping    |<-------- Scrub Stream (internal)
           |________________|
                   |
                   |
                   |
                   v
                 Product (Precipitation feed)

 NB. Extraction Feed (from decantation-filtration) goes to solvent extraction
 NB. Stripping Feed (from precipitation) goes to Stripping
 NB. Product Stream goes to Precipitation


   + Solvent Extraction
      0.1M Alamine 336 (TOA)
      Dilutent is Kerosene modified with 5% Isodecanol
      Aqueous/Organic = 3.0
      4 Mixer-Settler Units
      Settler Unit Area = 1400ft^2
      Total Setter Area = 5600ft^2
      
   + Scrubbing
      Utilizes Acidic Water as Wash Stream
   + Stripping
      Utilizes Acidic Sodium Chloride Solution

   Source of info:
   -https://www-pub.iaea.org/MTCD/Publications/PDF/trs359_web.pdf (pg. 189)
   -https://documents.deq.utah.gov/legacy/businesses/e/energy-fuels-resources-usa/docs/2007/05May/VOLUME%201.pdf (pg. 18)
"""

import logging

import math
from scipy.integrate import odeint
import numpy as np

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity
from cortix import Species

import unit

class Solvex(Module):
    """Solvent Extraction system.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: Decantation_Filtration, Precipitation.
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['extraction-feed', 'stripping-feed',
                                    'product', 'raffinate']

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
        '''
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
        '''

        # Initialization
        # Placeholder Values to Test Code
        
        self.extraction_feed_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.extraction_feed_mass_density = 1.0 * unit.kg / unit.liter

        self.extraction_raffinate_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.extraction_raffinate_mass_density = 1.0 * unit.kg / unit.liter

        self.extraction_product_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.extraction_product_mass_density = 1.0 * unit.kg / unit.liter

        self.stripping_feed_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.stripping_feed_mass_density = 1.0 * unit.kg / unit.liter

        self.stripping_product_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.stripping_product_mass_density = 1.0 * unit.kg / unit.liter

        # Derived quantities
        '''
        self.rho_p = 0.0
        self.rey_p = 0.0
        self.nusselt_p = 0.0

        self.nusselt_s = 0.0

        self.heat_sink_pwr = 0.0
        '''

        #***************************************************************************************
        # E X T R A C T I O N
        #***************************************************************************************

        # Extraction Feed Phase History (internal state/external)
        quantities = list()
        species = list()

        extraction_feed_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=self.extraction_feed_mass_flowrate,
                                          latex_name=r'$\dot{m}_1$',
                                          info='Extraction Feed Mass Flowrate')
        quantities.append(extraction_feed_mass_flowrate)

        extraction_feed_mass_density = Quantity(name='mass_density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=self.extraction_feed_mass_density,
                                         latex_name=r'$\rho$',
                                         info='Extraction Feed Mass Density')
        quantities.append(extraction_feed_mass_density)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(a)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O',formula_name='H2O(a)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        u6_aqu = Species( name='U-VI',formula_name='UO2^2+(a)',
                atoms=['U','2*O'],
                info='UO2$^{2+}$')
        species.append(u6_aqu)

        h2so4_feed = Species(name='H2-SO4',formula_name='H2SO4(a)',
                           atoms=['2*H','S','4*O'],
                           info='H2-SO4')
        species.append(h2so4_feed)

        hPlus_aqu = Species( name='H+',formula_name='H^+(a)',
                atoms=['H'],
                info='H$^+$')
        species.append(hPlus_aqu)
        
        toaso4 = Species(name='C24H51N-SO4',formula_name='C24H51NSO4(org)',
                         atoms=['24*C','51*H','N','S','4*O'],
                         info='C24H51N-SO4')
        species.append(toaso4)

        toauo2so43 = Species(name='C24H51N-UO2-(SO4)3',
                             formula_name='C24H51NUO2(SO4)3(org)',
                             atoms=['24*C','51*H','N','U','3*S','14*O'],
                             info='C24H51N-UO2-(SO4)3')
        species.append(toauo2so43)
        
        self.extraction_feed_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s', quantities=quantities, species=species)

        # Extraction Raffinate Phase History (internal state/external)
        quantities = list()
        species = list()

        extraction_raffinate_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=self.extraction_raffinate_mass_flowrate,
                                          latex_name=r'$\dot{m}_2$',
                                          info='Extraction Raffinate Mass Flowrate')
        quantities.append(extraction_raffinate_mass_flowrate)

        extraction_raffinate_mass_density = Quantity(name='mass_density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=self.extraction_raffinate_mass_density,
                                         latex_name=r'$\rho$',
                                         info='Extraction Raffinate Mass Density')
        quantities.append(extraction_raffinate_mass_density)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(a)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O',formula_name='H2O(a)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        u6_aqu = Species( name='U-VI',formula_name='UO2^2+(a)',
                atoms=['U','2*O'],
                info='UO2$^{2+}$')
        species.append(u6_aqu)

        h2so4_feed = Species(name='H2SO4',formula_name='H2SO4(a)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_feed)

        hPlus_aqu = Species( name='H+',formula_name='H^+(a)',
                atoms=['H'],
                info='H$^+$')
        species.append(hPlus_aqu)

        self.extraction_raffinate_phase = Phase(time_stamp=self.initial_time,
                                                time_unit='s', quantities=quantities, species=species)

        # Extraction Product Phase History (internal state/external)
        quantities = list()
        species = list()

        extraction_product_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=self.extraction_product_mass_flowrate,
                                          latex_name=r'$\dot{m}_3$',
                                          info='Extraction Product Mass Flowrate')
        quantities.append(extraction_product_mass_flowrate)

        extraction_product_mass_density = Quantity(name='mass_density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=self.extraction_product_mass_density,
                                         latex_name=r'$\rho$',
                                         info='Extraction Product Mass Density')
        quantities.append(extraction_product_mass_density)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(a)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O',formula_name='H2O(a)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        u6_aqu = Species( name='U-VI',formula_name='UO2^2+(a)',
                atoms=['U','2*O'],
                info='UO2$^{2+}$')
        species.append(u6_aqu)

        h2so4_feed = Species(name='H2-SO4',formula_name='H2SO4(a)',
                           atoms=['2*H','S','4*O'],
                           info='H2-SO4')
        species.append(h2so4_feed)

        hPlus_aqu = Species( name='H+',formula_name='H^+(a)',
                atoms=['H'],
                info='H$^+$')
        species.append(hPlus_aqu)
        
        toaso4 = Species(name='C24H51N-SO4',formula_name='C24H51NSO4(org)',
                         atoms=['24*C','51*H','N','S','4*O'],
                         info='C24H51N-SO4')
        species.append(toaso4)
 
        toauo2so43 = Species(name='C24H51N-UO2-(SO4)3',
                             formula_name='C24H51NUO2(SO4)3(org)',
                             atoms=['24*C','51*H','N','U','3*S','14*O'],
                             info='C24H51N-UO2-(SO4)3')
        species.append(toauo2so43)

        self.extraction_product_phase = Phase(time_stamp=self.initial_time,
                                              time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # S T R I P P I N G
        #***************************************************************************************

        # Stripping Feed Phase History (internal state/external)
        quantities = list()
        species = list()

        stripping_feed_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=self.stripping_feed_mass_flowrate,
                                          latex_name=r'$\dot{m}_4$',
                                          info='Stripping Feed Mass Flowrate')
        quantities.append(stripping_feed_mass_flowrate)

        stripping_feed_mass_density = Quantity(name='mass_density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=self.stripping_feed_mass_density,
                                         latex_name=r'$\rho$',
                                         info='Stripping Feed Mass Density')
        quantities.append(stripping_feed_mass_density)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(a)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O',formula_name='H2O(a)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        u6_aqu = Species( name='U-VI',formula_name='UO2^2+(a)',
                atoms=['U','2*O'],
                info='UO2$^{2+}$')
        species.append(u6_aqu)

        h2so4_feed = Species(name='H2-SO4',formula_name='H2SO4(a)',
                           atoms=['2*H','S','4*O'],
                           info='H2-SO4')
        species.append(h2so4_feed)

        hPlus_aqu = Species( name='H+',formula_name='H^+(a)',
                atoms=['H'],
                info='H$^+$')
        species.append(hPlus_aqu)
        
        toaso4 = Species(name='C24H51N-SO4',formula_name='C24H51NSO4(org)',
                         atoms=['24*C','51*H','N','S','4*O'],
                         info='C24H51N-SO4')
        species.append(toaso4)
        
        toauo2so43 = Species(name='C24H51N-UO2-(SO4)3',
                             formula_name='C24H51NUO2(SO4)3(org)',
                             atoms=['24*C','51*H','N','U','3*S','14*O'],
                             info='C24H51N-UO2-(SO4)3')
        species.append(toauo2so43)

        nh4oh = Species(name='NH4-OH',formula_name='NH4OH',
                        atoms=['N','5*H','O'],info='NH4-OH')
        species.append(nh4oh)
        
        toa = Species(name='C24H51N',formula_name='C24H51N',
                       atoms=['24*C','51*H','N'],info='C24H51N')
        species.append(toa)
        
        self.stripping_feed_phase = Phase(time_stamp=self.initial_time,
                                          time_unit='s', quantities=quantities, species=species)

        # Stripping Product Phase History (internal state/external)
        quantities = list()
        species = list()

        stripping_product_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=self.stripping_product_mass_flowrate,
                                          latex_name=r'$\dot{m}_5$',
                                          info='Stripping Product Mass Flowrate')
        quantities.append(stripping_product_mass_flowrate)

        stripping_product_mass_density = Quantity(name='mass_density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=self.stripping_product_mass_density,
                                         latex_name=r'$\rho$',
                                         info='Stripping Product Mass Density')
        quantities.append(stripping_product_mass_density)

        uo2so434minus_feed = Species(name='UO2-(SO4)3^4-',formula_name='UO2(SO4)3^4-(a)',
                           atoms=['U','2*O','3*S','12*O'],
                           info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_feed)

        h2o_feed = Species(name='H2O',formula_name='H2O(a)',
                           atoms=['2*H','O'],
                           info='H2O')
        species.append(h2o_feed)

        u6_aqu = Species( name='U-VI',formula_name='UO2^2+(a)',
                atoms=['U','2*O'],
                info='UO2$^{2+}$')
        species.append(u6_aqu)

        h2so4_feed = Species(name='H2SO4',formula_name='H2SO4(a)',
                           atoms=['2*H','S','4*O'],
                           info='H2SO4')
        species.append(h2so4_feed)

        hPlus_aqu = Species( name='H+',formula_name='H^+(a)',
                atoms=['H'],
                info='H$^+$')
        species.append(hPlus_aqu)

        self.stripping_product_phase = Phase(time_stamp=self.initial_time,
                                             time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # S T A T E  P H A S E
        #***************************************************************************************

        '''
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
        '''

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

        self.end_time = time # correct the final time if needed

    def __call_ports(self, time):

        # Interactions in the uranium-inflow port
        #----------------------------------------
        # One way "from" extraction-feed

        # Receive from
        if self.get_port('extraction-feed').connected_port:

            self.send(time, 'extraction-feed')

            (check_time, extraction_feed) = self.recv('extraction-feed')
            assert abs(check_time-time) <= 1e-6

            '''
            self.primary_inflow_temp = primary_inflow['temperature']
            self.primary_ressure = primary_inflow['pressure']
            self.primary_mass_flowrate = primary_inflow['mass_flowrate']
            '''

        # Interactions in the secondary-inflow port
        #----------------------------------------
        # One way "from" stripping-feed

        # Receive from
        if self.get_port('stripping-feed').connected_port:

            self.send(time, 'stripping-feed')

            (check_time, stripping_feed) = self.recv('striping-feed')
            assert abs(check_time-time) <= 1e-6

            '''
            self.secondary_inflow_temp = secondary_inflow['temperature']
            self.secondary_pressure = secondary_inflow['pressure']
            self.secondary_mass_flowrate = secondary_inflow['mass_flowrate']
            '''

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # One way "to" product

        # Send to
        if self.get_port('product').connected_port:

            msg_time = self.recv('product')

            product = dict()
            product['mass-flowrate'] = self.stripping_product_phase.get_value('mass-flowrate',msg_time)
            product['mass-density'] = self.stripping_product_phase.get_value('mass-density',msg_time)
            '''
            product['temperature'] = temp
            product['pressure'] = self.primary_pressure
            product['mass_flowrate'] = self.primary_mass_flowrate
            product['quality'] = 0.0
            '''

            self.send((msg_time, product), 'product')

        # Interactions in the secondary-outflow port
        #-----------------------------------------
        # One way "to" raffinate

        # Send to
        if self.get_port('raffinate').connected_port:

            msg_time = self.recv('raffinate')

            raffinate = dict()
            product['mass-flowrate'] = self.extraction_raffinate_phase.get_value('mass-flowrate',msg_time)
            product['mass-density'] = self.extraction_raffinate_phase.get_value('mass-density',msg_time)
            '''
            raffinate['temperature'] = temp
            raffinate['pressure'] = press
            raffinate['mass_flowrate'] = flowrate
            raffinate['total_heat_power'] = -self.heat_sink_pwr
            '''

            self.send((msg_time, raffinate), 'raffinate')

    def __step(self, time=0.0):
        """Stepping Solvex in time
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

        #Time Step with Constant Value to Test Code
        tmp = self.extraction_feed_phase.get_row(time)
        print(tmp)
        mass_flowrate = self.extraction_feed_phase.get_value('mass_flowrate', time)

        time += self.time_step
        
        self.extraction_feed_phase.add_row(time, tmp)

        self.extraction_feed_phase.set_value('mass_flowrate', mass_flowrate, time)
        #self.extraction_feed_phase.set_value('mass_density', rho_preleach, time)


        '''
        self.primary_outflow_phase.add_row(time, primary_outflow)
        self.primary_outflow_phase.set_value('temp', temp_p, time)
        self.primary_outflow_phase.set_value('flowrate', self.primary_mass_flowrate, time)

        self.secondary_inflow_phase.add_row(time, secondary_inflow)
        self.secondary_inflow_phase.set_value('temp', self.secondary_inflow_temp, time)
        self.secondary_inflow_phase.set_value('flowrate', self.secondary_mass_flowrate, time)

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
        '''

        return time
