#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module
   This module is a model of the Solvent Extraction process in the White Mesa Uranium Milling Plant


                              |
                              |  Extraction Feed (from Decantation-Filtration filtrate)
                              |
                              v
                      |----------------|
                      |                |
 Organic Feed ------->|                |------> Organic Product (to Scrubbing internal)
 (internal source)    |    Solvent     |
                      |   Extraction   |
 Raffinate <----------|                |<------ Scrub Raffinate (from Scrubbing internal)
 Stream (to CCD Bank) |                |
                      |----------------|
  Organic Product <---|                |<------ Organic Feed (Organic Product from Solv. Extr. internal)
  (to Strip internal) |                |
                      |   Scrubbing    |
   Scrub Feed ------->|                |-------> Scrub Raffinate (to Solvent Extraction internal)
   (internal source)  |                |
                      |----------------|
  Organic Feed ------>|                |-------> Organic Regeneration (to Solvent Extraction not done)
(from Scrub internal) |   Stripping    |
                      |                |<------ Stripping Feed (internal source)
                      |________________|<------ Stripping Feed (from Precipitation not implemented)
                              |
                              |
                              |
                              v
                       Stripping Product (Precipitation feed)

 NB. Extraction Feed (from decantation-filtration) goes to solvent extraction

   + Solvent Extraction
      - 4 Mixer-settlers units. Area per unit 1400 ft^2, height 6 ft
      - Flowrates:
      - Aqueous/Organic = 3.0
      - 0.1M Alamine 336 (tri-octylamine TOA)
      - Dilutent is kerosene modified with 5% isodecanol
      - Solvent mixture: 2.5% TOA, 2.5% isodecanol, and 95% kerosene

   + Scrubbing
      - 1 Mixer-settlers units. Area per unit 1400 ft^2, height 6 ft

   + Stripping
      - 4 Mixer-settlers units. Area per unit 1400 ft^2, height 6 ft

   + Solvent Regeneration
      - 1 Mixer-settlers units. Area per unit 1400 ft^2, height 6 ft

   Source of info:
   -https://www-pub.iaea.org/MTCD/Publications/PDF/trs359_web.pdf (pg. 189, 337)
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
from cortix import Units as unit

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

        # Extraction
        self.extract_tank_volume = 4 * 1400 * unit.ft**2 * 6 * unit.ft
        self.aqueous_to_organic_volume_ratio = 3
        self.aqueous_raffinate_to_feed_rho_ratio = 0.73
        self.organic_product_to_feed_rho_ratio = 1.25

        # Scrubbing
        self.scrub_tank_volume = 1 * 1400 * unit.ft**2 * 6 * unit.ft

        # Stripping
        self.strip_tank_volume = 4 * 1400 * unit.ft**2 * 6 * unit.ft

        # Initialization

        # Solvent Extraction
        if self.get_port('extraction-feed').connected_port:
            self.extract_feed_mass_flowrate = 0.0 * unit.liter / unit.minute
            self.extract_feed_mass_density = 0.0 * unit.kg / unit.liter
        else:
            self.extract_feed_mass_flowrate = 100.0 * unit.kg / unit.minute
            self.extract_feed_mass_density = 1.6 * unit.kg / unit.liter

        self.extract_organic_feed_mass_flowrate = self.extract_feed_mass_flowrate /\
                                                 self.aqueous_to_organic_volume_ratio

        self.extract_organic_feed_mass_density = 0.8 * unit.kg/unit.liter

        #***************************************************************************************
        # E X T R A C T I O N
        #***************************************************************************************

        # Extraction Raffinate Phase History (internal state/external)
        quantities = list()
        species = list()

        extraction_raffinate_mass_flowrate = Quantity(name='mass-flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=0.0,
                                          latex_name=r'$\dot{m}_{r}$',
                                          info='Extraction Raffinate Mass Flowrate')
        quantities.append(extraction_raffinate_mass_flowrate)

        extraction_raffinate_mass_density = Quantity(name='mass-density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=0.0,
                                         latex_name=r'$\rho_{r}$',
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

        self.extract_raffinate_phase = Phase(time_stamp=self.initial_time,
                                            time_unit='s', quantities=quantities, species=species)

        # Extraction Product Phase History (internal state/external)
        quantities = list()
        species = list()

        extraction_product_mass_flowrate = Quantity(name='mass-flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=0.0,
                                          latex_name=r'$\dot{m}_3$',
                                          info='Extraction Product Mass Flowrate')
        quantities.append(extraction_product_mass_flowrate)

        extraction_product_mass_density = Quantity(name='mass-density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=0.0,
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

        self.extract_product_phase = Phase(time_stamp=self.initial_time,
                                          time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # S C R U B B I N G
        #***************************************************************************************

        # Scrub Raffinate Phase History (internal state/external)
        quantities = list()
        species = list()

        extraction_raffinate_mass_flowrate = Quantity(name='mass-flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=0.0,
                                          latex_name=r'$\dot{m}_2$',
                                          info='Scrubbing Raffinate Mass Flowrate')
        quantities.append(extraction_raffinate_mass_flowrate)

        extraction_raffinate_mass_density = Quantity(name='mass-density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=0.0,
                                         latex_name=r'$\rho$',
                                         info='Scrubbing Raffinate Mass Density')
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

        self.scrub_raffinate_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s', quantities=quantities, species=species)

        #***************************************************************************************
        # S T R I P P I N G
        #***************************************************************************************

        # Stripping Feed Phase History (internal state/external)
        quantities = list()
        species = list()

        stripping_feed_mass_flowrate = Quantity(name='mass-flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=0.0,
                                          latex_name=r'$\dot{m}_4$',
                                          info='Stripping Feed Mass Flowrate')
        quantities.append(stripping_feed_mass_flowrate)

        stripping_feed_mass_density = Quantity(name='mass-density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=0.0,
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

        stripping_product_mass_flowrate = Quantity(name='mass-flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=0.0,
                                          latex_name=r'$\dot{m}_5$',
                                          info='Stripping Product Mass Flowrate')
        quantities.append(stripping_product_mass_flowrate)

        stripping_product_mass_density = Quantity(name='mass-density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=0.0,
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

        # Solvent Extraction State
        quantities = list()
        species = list()

        liq_volume = Quantity(name='aqueous-volume',
                        formal_name='v', unit='m$^3$',
                        value=0.0,
                        latex_name=r'$V_{a}$',
                        info='Solvent Extraction Aqueous Volume')
        quantities.append(liq_volume)

        liq_volume = Quantity(name='organic-volume',
                        formal_name='v', unit='m$^3$',
                        value=0.0,
                        latex_name=r'$V_{o}$',
                        info='Solvent Extraction Organic Volume')
        quantities.append(liq_volume)

        liq_volume = Quantity(name='liquid-volume',
                        formal_name='v', unit='m$^3$',
                        value=0.0,
                        latex_name=r'$V_{oa}$',
                        info='Solvent Extraction Liquid (Aqueous + Organic)  Volume')
        quantities.append(liq_volume)

        self.extract_state_phase = Phase(time_stamp=self.initial_time, time_unit='s',
                                        quantities=quantities, species=species)

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

        # Interactions in the uranium-inflow port
        #----------------------------------------
        # One way "from" extraction-feed

        # Receive from
        if self.get_port('extraction-feed').connected_port:

            self.send(time, 'extraction-feed')

            (check_time, extraction_feed) = self.recv('extraction-feed')
            assert abs(check_time-time) <= 1e-6

        # Interactions in the secondary-inflow port
        #----------------------------------------
        # One way "from" stripping-feed

        # Receive from
        if self.get_port('stripping-feed').connected_port:

            self.send(time, 'stripping-feed')

            (check_time, stripping_feed) = self.recv('striping-feed')
            assert abs(check_time-time) <= 1e-6

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # One way "to" product

        # Send to
        if self.get_port('product').connected_port:

            msg_time = 0.0
            #msg_time = self.recv('product')
            #if msg_time is None:
            #    msg_time = 0.0

            product = dict()
            product['mass-flowrate'] = self.stripping_product_phase.get_value('mass-flowrate',msg_time)
            product['mass-density'] = self.stripping_product_phase.get_value('mass-density',msg_time)

            self.send((msg_time, product), 'product')

        # Interactions in the secondary-outflow port
        #-----------------------------------------
        # One way "to" raffinate

        # Send to
        if self.get_port('raffinate').connected_port:

            msg_time = self.recv('raffinate')

            raffinate = dict()
            raffinate['mass-flowrate'] = self.extraction_raffinate_phase.get_value('mass-flowrate',msg_time)
            raffinate['mass-density'] = self.extraction_raffinate_phase.get_value('mass-density',msg_time)

            self.send((msg_time, raffinate), 'raffinate')

    def __step(self, time=0.0):
        """Stepping Extraction in time
        """

        #------------------------------------
        # Evolve the Solvent Extraction State
        #------------------------------------
        # Ideal flow mixing

        # Aqueous
        extract_raffinate_mass_flowrate_initial = self.extract_raffinate_phase.get_value('mass-flowrate', time)

        scrub_raffinate_mass_flowrate_inflow = self.scrub_raffinate_phase.get_value('mass-flowrate', time)
        scrub_raffinate_mass_density_inflow = self.scrub_raffinate_phase.get_value('mass-density', time)

        aq_feed_mass_flowrate = self.extract_feed_mass_flowrate
        aq_rho_feed = self.extract_feed_mass_density

        if aq_feed_mass_flowrate == 0.0:
            aq_feed_vol_flowrate = 0.0
        else:
            aq_feed_vol_flowrate = aq_feed_mass_flowrate/aq_rho_feed

        # Ideal solution
        aqueous_mass_flowrate_inflow = aq_feed_mass_flowrate + scrub_raffinate_mass_flowrate_inflow
        aqueous_rho = (aq_rho_feed * self.aqueous_raffinate_to_feed_rho_ratio + \
                      scrub_raffinate_mass_density_inflow) / 2
        aqueous_vol_flowrate_inflow = aqueous_mass_flowrate_inflow/aqueous_rho

        aqueous_mass_flowrate_initial = extract_raffinate_mass_flowrate_initial
        aqueous_vol_flowrate_initial = aqueous_mass_flowrate_initial/aqueous_rho

        aqueous_volume_initial = self.extract_state_phase.get_value('aqueous-volume', time)

        aqueous_volume = aqueous_volume_initial + \
                         (aqueous_vol_flowrate_inflow - aqueous_vol_flowrate_initial) * self.time_step

        # Organic
        extract_organic_mass_flowrate_initial = self.extract_product_phase.get_value('mass-flowrate', time)

        org_feed_mass_flowrate = self.extract_organic_feed_mass_flowrate

        org_rho_feed = self.extract_organic_feed_mass_density

        #if org_feed_mass_flowrate == 0.0:
        #    org_feed_vol_flowrate = 0.0
        #else:
        #    org_feed_vol_flowrate = org_feed_mass_flowrate/org_rho_feed

        # Ideal solution
        organic_mass_flowrate_inflow = org_feed_mass_flowrate
        organic_rho = org_rho_feed * self.organic_product_to_feed_rho_ratio
        organic_vol_flowrate_inflow = organic_mass_flowrate_inflow/organic_rho

        organic_mass_flowrate_initial = extract_organic_mass_flowrate_initial
        organic_vol_flowrate_initial = organic_mass_flowrate_initial/organic_rho

        organic_volume_initial = self.extract_state_phase.get_value('organic-volume', time)

        organic_volume = organic_volume_initial + \
                         (organic_vol_flowrate_inflow - organic_vol_flowrate_initial) * self.time_step

        liquid_volume = aqueous_volume + organic_volume

        # Place-holder for mass balance
        if liquid_volume >= 0.5 * self.extract_tank_volume:

            flow_residence_time = self.extract_tank_volume /\
                                  (aqueous_vol_flowrate_inflow + organic_vol_flowrate_inflow)

            mass_flowrate_raffinate = aqueous_mass_flowrate_inflow + \
                                      math.exp(-self.time_step/flow_residence_time) * \
                                      (aqueous_mass_flowrate_initial - aqueous_mass_flowrate_inflow)

            mass_flowrate_organic = organic_mass_flowrate_inflow + \
                                    math.exp(-self.time_step/flow_residence_time) * \
                                    (organic_mass_flowrate_initial - organic_mass_flowrate_inflow)
        else:
            mass_flowrate_raffinate = 0.0
            mass_flowrate_organic   = 0.0


        tmp_extract_state = self.extract_state_phase.get_row(time)
        tmp_extract_raffinate = self.extract_raffinate_phase.get_row(time)
        tmp_extract_product = self.extract_product_phase.get_row(time)

        #---------------------------
        # Evolve the Scrubbing State
        #---------------------------
        # Ideal flow mixing

        tmp_scrub_raffinate = self.extract_raffinate_phase.get_row(time)

        #----------------------------
        # Step All Quantities in Time
        #----------------------------

        time += self.time_step

        #-------
        # SolvEx
        #-------
        self.extract_state_phase.add_row(time, tmp_extract_state)
        self.extract_state_phase.set_value('aqueous-volume', aqueous_volume, time)
        self.extract_state_phase.set_value('organic-volume', organic_volume, time)
        self.extract_state_phase.set_value('liquid-volume', liquid_volume, time)

        self.extract_raffinate_phase.add_row(time, tmp_extract_raffinate)
        self.extract_raffinate_phase.set_value('mass-flowrate', mass_flowrate_raffinate, time)
        self.extract_raffinate_phase.set_value('mass-density', aqueous_rho, time)

        self.extract_product_phase.add_row(time, tmp_extract_product)
        self.extract_product_phase.set_value('mass-flowrate', mass_flowrate_organic, time)
        self.extract_product_phase.set_value('mass-density', organic_rho, time)


        self.scrub_raffinate_phase.add_row(time, tmp_scrub_raffinate)
        self.scrub_raffinate_phase.set_value('mass-flowrate', scrub_raffinate_mass_flowrate_inflow)

        return time
