#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix Module.
   Leaching process in the White Mesa Milling Plant.


                  Wet Ore
                 (internal)
                     |
                     |
                     |
                     v
             |----------------|
             |                |
             |  Pre-leaching  |<-------- Pre-Leach Feed (CCD overflow from Decantation Module)
             |                |
             |                |
             |                |<-------- STD Underflow (from single-tank Decantation Module)
             |  Acid-leaching |
             |                |<-------- Acids (internal) H2S04, NaCI03, Steam)
             |________________|
                  |       |
                  |       |
                  |       |
                  |       v
                  |    Pre-leach Product (STD Single-Tank Decantation)
                  v
    Acid-leach Product (CCD Decantation)


   + Pre-Leaching

   Add info here... what ore mineral (brannerite)?
                    what oxidation process?

                    Carnotite sandstone 0.2% U3O8, 1.5-2.0% V2O5
                    Arizona Strip breccia pipe 0.5-0.9% U3O8
                    Example of Breccia Pipe Ore Metal Concentrations: 3000 ppm Arsenic, 200 ppm Cobalt,
                    8000 ppm copper, 6 ppm mercury, 260 ppm molybdenum, 500 ppm nickel, 1% lead,
                    3000 ppm uranium (0.3%), 150 ppm Zinc. These concentrations can vary from 0.2-2% (20000-2000 ppm)
                    The original mill design planned for 0.2-0.9% uranium. This is a relatively high concentration
                    compared to many mines but is well within the averages/usuals for most mines.
                    Uranium typically exists in the ores in the form of U3O8.

      *Pre-Leach Ore Feed
          -Mix of 1ton of ore and water
          -55-58% solids. On a basis of 1 ton of ore feed....
          -1000 kg ore
          -1000 / 0.55 = 1818.18 kg total
          -818.18 kg water

      *Pre-Leach Output
          -22% solids
          -
      Pulp or solids density is 22%. So 1t of ore into the preleach leaves with 220kg of solids.
   + Acid-Leaching
      *Chemistry EQNS
          -Uranium
              Typical oxidation of uranium from solid ore
              1.  UO3(s)+2H^+(aq) --> UO2^2+(aq) + H2O(aq)
          -Uranium and Iron
              2.  UO2(s) +2Fe^3+(sq) --> UO2^2+(aq) + 2Fe^2+(aq)
          -Uranium and Sulfiric Acid
              3.  UO2^2+(aq) + 3(SO4^2-)(aq) --> (UO2(SO4)3)^4-(aq)
              4.  OR UO2^2+(aq) + 3SO4(2-)(aq) --> UO2(SO4)2^2-(aq)
          -Gold/Au
                                        Inert
          -Sodium Chloride (Dissolution)
              5.  NaClO3(aq) --> ClO3^-(aq) + Na+(aq)
          -Iron and Sulfiric Acid
              6.  Fe3O4 + H2SO4 --> FeSO4 + Fe2O3 +H2O
              7.  2Fe^2+(aq) + 1/3ClO3^2-(aq) + 2H^+(aq) --> 2Fe^3+(aq) + 1/3Cl^-(aq) + H2O(aq)
          -Copper
      - Capacity: 1 t of ore
      - Acid (H2SO4) amount: 20 kg/t ore
          Optimal concentration of the acid feed for selectivity seems to be 10-48% sulfuric acid for the aqueous solution
      - Temperature: 40C
          Research paper showed that 40C was preferred to 50C, 60C, and 80C for selectivity purposes. Lower temp might be better
      - Residual H2SO4: 50 g/L free acid

   Source of info:
      https://pubs.usgs.gov/sir/2010/5025/pdf/sir2010-5025_availability.pdf
      https://www-pub.iaea.org/MTCD/publications/PDF/TE_1629_web.pdf
      - Used to find chemical equations for Uranium+Look at kinetics
      https://www.sciencedirect.com/science/article/pii/S1738573321005970
      - Kinetic Equations?
      https://repository.up.ac.za/bitstream/handle/2263/61336/Sililo_Modelling_2017.pdf?sequence=1
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

class Leaching(Module):
    """Heap Leach.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: Filtration/Decantation.
    See instance attribute `port_names_expected`.
    """

    def __init__(self):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['pre-leach-feed', 'acid-leach-feed',
                                    'pre-leach-product', 'acid-leach-product']

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

        self.preleach_vol = 45 * unit.meter**3
        self.wet_ore_feed_mass_density = 5000 * unit.kg/unit.meter**3

        # Initialization

        # Pre-leaching [These values are temporary, real ones will have to be added]
        self.wet_ore_feed_mass_flowrate = 2727 * unit.kg/unit.minute
        self.wet_ore_feed_solid_mass_fraction = 55/100

        self.wet_ore_mass_density = 1.0 * unit.kg / unit.liter
        self.wet_ore_solids_massfrac = 100 * unit.ppm

        self.preleach_feed_mass_flowrate = 3500 * unit.kg / unit.minute
        self.preleach_feed_mass_density = 1.6 * unit.kg / unit.liter
        self.preleach_feed_solids_massfrac = 100 * unit.ppm

        self.preleach_output_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.preleach_output_mass_density = 1.0 * unit.kg / unit.liter
        self.preleach_output_solids_massfrac = 100 * unit.ppm

        # Acid-leaching
        self.std_underflow_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.std_underflow_mass_density = 1.0 * unit.kg / unit.liter
        self.std_underflow_solids_massfrac = 100 * unit.ppm

        self.acids_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.acids_mass_density = 1.0 * unit.kg / unit.liter

        self.acid_leach_output_mass_flowrate = 1.0 * unit.liter / unit.minute
        self.acid_leach_output_mass_density = 1.0 * unit.kg / unit.liter
        self.acid_leach_output_solids_massfrac = 100 * unit.ppm

        # ***************************************************************************************
        # P R E - L E A C H I N G
        # ***************************************************************************************

        # Wet Ore Feed Phase History (internal)
        quantities = list()
        species = list()

        wet_ore_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=self.preleach_feed_mass_flowrate,
                                          latex_name=r'$\dot{m}_1$',
                                          info='Wet Ore Mass Flowrate')
        quantities.append(wet_ore_mass_flowrate)

        wet_ore_mass_density = Quantity(name='mass_density',
                                         formal_name='rho', unit='kg/m^3',
                                         value=self.wet_ore_mass_density,
                                         latex_name=r'$\rho$',
                                         info='Wet Ore Feed Mass Density')
        quantities.append(wet_ore_mass_density)

        wet_ore_solids_massfrac = Quantity(name='solids_massfrac',
                                            formal_name='solids_massfrac', unit='ppm',
                                            value=self.wet_ore_solids_massfrac,
                                            latex_name=r'$C_1$',
                                            info='Wet Ore Solids Mass Fraction')

        quantities.append(wet_ore_solids_massfrac)

        uo2so434minus_wet_ore = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(aq)',
                                         atoms=['U', '2*O', '3*S', '12*O'],
                                         info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_wet_ore)

        h2o_wet_ore = Species(name='H2O', formula_name='H2O(aq)',
                               atoms=['2*H', 'O'],
                               info='H2O')
        species.append(h2o_wet_ore)

        h2so4_wet_ore = Species(name='H2SO4', formula_name='H2SO4(aq)',
                                 atoms=['2*H', 'S', '4*O'],
                                 info='H2SO4')
        species.append(h2so4_wet_ore)

        iron_wet_ore = Species(name='Fe', formula_name='Fe(s)',
                                atoms=['Fe'],
                                info='Fe')
        species.append(iron_wet_ore)

        copper_wet_ore = Species(name='Cu', formula_name='Cu(s)',
                                  atoms=['Cu'],
                                  info='Cu')
        species.append(copper_wet_ore)

        gold_wet_ore = Species(name='Au', formula_name='Au(s)',
                                atoms=['Au'],
                                info='Au')
        species.append(gold_wet_ore)

        self.wet_ore_feed_phase = Phase(time_stamp=self.initial_time,
                                        time_unit='s', quantities=quantities, species=species)

        # Pre-Leach Phase History (Goes to Single-Tank Decantation)
        quantities = list()
        species = list()

        preleach_output_mass_flowrate = Quantity(name='mass_flowrate',
                                          formal_name='mdot', unit='kg/s',
                                          value=0.0,
                                          latex_name=r'$\dot{m}_{p}$',
                                          info='Pre-Leach Mass Flowrate')
        quantities.append(preleach_output_mass_flowrate)

        preleach_output_mass_density = Quantity(name='mass_density',
                                     formal_name='rho', unit='kg/m^3',
                                     value=0.0,
                                     latex_name=r'$\rho$',
                                     info='Pre-Leach Mass Density')
        quantities.append(preleach_output_mass_density)

        preleach_output_solids_massfrac = Quantity(name='solids_massfrac',
                                            formal_name='solids_massfrac', unit='ppm',
                                            value=self.preleach_output_solids_massfrac,
                                            latex_name=r'$C_1$',
                                            info='Preleach Solids Mass Fraction')

        quantities.append(preleach_output_solids_massfrac)

        uo2so434minus_preleach_output = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(aq)',
                                         atoms=['U', '2*O', '3*S', '12*O'],
                                         info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_preleach_output)

        h2o_preleach_output = Species(name='H2O', formula_name='H2O(aq)',
                               atoms=['2*H', 'O'],
                               info='H2O')
        species.append(h2o_preleach_output)

        h2so4_preleach_output = Species(name='H2SO4', formula_name='H2SO4(aq)',
                                 atoms=['2*H', 'S', '4*O'],
                                 info='H2SO4')
        species.append(h2so4_preleach_output)

        iron_preleach_output = Species(name='Fe', formula_name='Fe(s)',
                                atoms=['Fe'],
                                info='Fe')
        species.append(iron_preleach_output)

        copper_preleach_output = Species(name='Cu', formula_name='Cu(s)',
                                  atoms=['Cu'],
                                  info='Cu')
        species.append(copper_preleach_output)

        gold_preleach_output = Species(name='Au', formula_name='Au(s)',
                                atoms=['Au'],
                                info='Au')
        species.append(gold_preleach_output)

        self.preleach_phase = Phase(time_stamp=self.initial_time,
                                            time_unit='s', quantities=quantities, species=species)

        # ***************************************************************************************
        # A C I D - L E A C H I N G
        # ***************************************************************************************

        # Acid-Leach Feed Phase History (STD underflow)
        quantities = list()
        species = list()

        underflow_mass_flowrate = Quantity(name='mass_flowrate',
                                           formal_name='mdot', unit='kg/s',
                                           value=self.std_underflow_mass_flowrate,
                                           latex_name=r'$\dot{m}_4$',
                                           info='STD Underflow Mass Flowrate')
        quantities.append(underflow_mass_flowrate)

        std_underflow_mass_density = Quantity(name='mass_density',
                                     formal_name='rho', unit='kg/m^3',
                                     value=self.std_underflow_mass_density,
                                     latex_name=r'$\rho$',
                                     info='STD Underflow Feed Mass Density')
        quantities.append(std_underflow_mass_density)

        underflow_solids_massfrac = Quantity(name='solids_massfrac',
                                             formal_name='solids_massfrac', unit='ppm',
                                             value=self.std_underflow_solids_massfrac,
                                             latex_name=r'$C_1$',
                                             info='STD Underflow Solids Mass Fraction')
        quantities.append(underflow_solids_massfrac)

        uo2so434minus_underflow = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(aq)',
                                          atoms=['U', '2*O', '3*S', '12*O'],
                                          info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_underflow)

        h2o_underflow = Species(name='H2O', formula_name='H2O(aq)',
                                atoms=['2*H', 'O'],
                                info='H2O')
        species.append(h2o_underflow)

        h2so4_underflow = Species(name='H2SO4', formula_name='H2SO4(aq)',
                                  atoms=['2*H', 'S', '4*O'],
                                  info='H2SO4')
        species.append(h2so4_underflow)

        iron_underflow = Species(name='Fe', formula_name='Fe(s)',
                                 atoms=['Fe'],
                                 info='Fe')
        species.append(iron_underflow)

        copper_underflow = Species(name='Cu', formula_name='Cu(s)',
                                   atoms=['Cu'],
                                   info='Cu')
        species.append(copper_underflow)

        gold_underflow = Species(name='Au', formula_name='Au(s)',
                                 atoms=['Au'],
                                 info='Au')
        species.append(gold_underflow)

        self.acidleach_feed_phase = Phase(time_stamp=self.initial_time,
                                          time_unit='s', quantities=quantities, species=species)

        # Acids feed phase history
        quantities = list()
        species = list()

        acids_mass_flowrate = Quantity(name='mass_flowrate',
                                           formal_name='mdot', unit='kg/s',
                                           value=self.acids_mass_flowrate,
                                           latex_name=r'$\dot{m}_4$',
                                           info='Acid Feed Mass Flowrate')
        quantities.append(acids_mass_flowrate)

        acids_mass_density = Quantity(name='mass_density',
                                              formal_name='rho', unit='kg/m^3',
                                              value=self.acids_mass_density,
                                              latex_name=r'$\rho$',
                                              info='Acids Feed Mass Density')
        quantities.append(acids_mass_density)

        nacio3_acids = Species(name='NaCIO3', formula_name='NaCIO3(aq)',
                                          atoms=['Na', 'C', 'I', '3*O'],
                                          info='NaCIO3')
        species.append(nacio3_acids)

        h2o_acids = Species(name='H2O', formula_name='H2O(g)',
                                atoms=['2*H', 'O'],
                                info='H2O')
        species.append(h2o_acids)

        h2so4_acids = Species(name='H2SO4', formula_name='H2SO4(aq)',
                                  atoms=['2*H', 'S', '4*O'],
                                  info='H2SO4')
        species.append(h2so4_acids)

        self.acids_phase = Phase(time_stamp=self.initial_time,
                                         time_unit='s', quantities=quantities, species=species)

        # Acid-Leach Product Phase History (Decantation feed)
        quantities = list()
        species = list()

        acid_leach_output_mass_flowrate = Quantity(name='mass_flowrate',
                                           formal_name='mdot', unit='kg/s',
                                           value=self.acid_leach_output_mass_flowrate,
                                           latex_name=r'$\dot{m}_4$',
                                           info='Acid Leach Output Mass Flowrate')
        quantities.append(acid_leach_output_mass_flowrate)

        acid_leach_output_mass_density = Quantity(name='mass_density',
                                              formal_name='rho', unit='kg/m^3',
                                              value=self.acid_leach_output_mass_density,
                                              latex_name=r'$\rho$',
                                              info='Acid Leach Output Mass Density')
        quantities.append(acid_leach_output_mass_density)

        acid_leach_output_solids_massfrac = Quantity(name='solids_massfrac',
                                             formal_name='solids_massfrac', unit='ppm',
                                             value=self.acid_leach_output_solids_massfrac,
                                             latex_name=r'$C_1$',
                                             info='Acid Leach Output Solids Mass Fraction')
        quantities.append(acid_leach_output_solids_massfrac)

        uo2so434minus_acid_leach_output = Species(name='UO2-(SO4)3^4-', formula_name='UO2(SO4)3^4-(aq)',
                                          atoms=['U', '2*O', '3*S', '12*O'],
                                          info='UO2-(SO4)3^4-')
        species.append(uo2so434minus_acid_leach_output)

        h2o_acid_leach_output = Species(name='H2O', formula_name='H2O(aq)',
                                atoms=['2*H', 'O'],
                                info='H2O')
        species.append(h2o_acid_leach_output)

        h2so4_acid_leach_output = Species(name='H2SO4', formula_name='H2SO4(aq)',
                                  atoms=['2*H', 'S', '4*O'],
                                  info='H2SO4')
        species.append(h2so4_acid_leach_output)

        iron_acid_leach_output = Species(name='Fe', formula_name='Fe(s)',
                                 atoms=['Fe'],
                                 info='Fe')
        species.append(iron_acid_leach_output)

        copper_acid_leach_output = Species(name='Cu', formula_name='Cu(s)',
                                   atoms=['Cu'],
                                   info='Cu')
        species.append(copper_acid_leach_output)

        gold_acid_leach_output = Species(name='Au', formula_name='Au(s)',
                                 atoms=['Au'],
                                 info='Au')
        species.append(gold_acid_leach_output)

        self.acidleach_product_phase = Phase(time_stamp=self.initial_time,
                                             time_unit='s', quantities=quantities, species=species)

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

        # Interactions in the pre-leach-feed port
        # Receive from
        if self.get_port('pre-leach-feed').connected_port:

            self.send(time, 'pre-leach-feed')

            (check_time, preleach_feed) = self.recv('pre-leach-feed')
            assert abs(check_time-time) <= 1e-6

            self.preleach_feed_mass_flowrate = preleach_feed['mass-flowrate']
            '''
            self.primary_ressure = primary_inflow['pressure']
            self.primary_mass_flowrate = primary_inflow['mass_flowrate']
            '''

        # Interactions in the acid-leach-feed port
        #----------------------------------------
        # One way "from" acid-leach-feed port

        # Receive from
        if self.get_port('acid-leach-feed').connected_port:

            self.send(time, 'acid-leach-feed')

            (check_time, stripping_feed) = self.recv('acid-leach-feed')
            assert abs(check_time-time) <= 1e-6

            '''
            self.secondary_inflow_temp = secondary_inflow['temperature']
            self.secondary_pressure = secondary_inflow['pressure']
            self.secondary_mass_flowrate = secondary_inflow['mass_flowrate']
            '''

        # Interactions in the pre-leach-product port
        #-----------------------------------------
        # One way "to" pre-leach-product port

        # Send to
        if self.get_port('pre-leach-product').connected_port:

            msg_time = self.recv('pre-leach-product')

            product = dict()
            product['mass-flowrate'] = self.preleach_phase.get_value('mass-flowrate',msg_time)
            product['mass-density'] = self.preleach_phase.get_value('mass-density',msg_time)
            '''
            product['temperature'] = temp
            product['pressure'] = self.primary_pressure
            product['quality'] = 0.0
            '''

            self.send((msg_time, product), 'pre-leach-product')

        # Interactions in the acid-leach-product port
        #-----------------------------------------
        # One way "to" acid-leach-product-port

        # Send to
        if self.get_port('acid-leach-product').connected_port:

            msg_time = self.recv('acid-leach-product')

            #temp = self.xxxx.get_value('temp', msg_time)

            raffinate = dict()
            '''
            raffinate['temperature'] = temp
            raffinate['pressure'] = press
            raffinate['mass_flowrate'] = flowrate
            raffinate['total_heat_power'] = -self.heat_sink_pwr
            '''

            self.send((msg_time, raffinate), 'acid-leach-product')

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

        # Evolving the pre-leach state
        mass_flowrate_initial = self.preleach_phase.get_value('mass_flowrate', time)

        wet_ore_mass_flowrate = self.wet_ore_feed_mass_flowrate
        preleach_feed_mass_flowrate = self.preleach_feed_mass_flowrate
        mass_flowrate_inflow = wet_ore_mass_flowrate + preleach_feed_mass_flowrate

        rho_preleach_feed = self.preleach_feed_mass_density
        rho_wet_ore = self.wet_ore_mass_density

        rho_preleach = rho_preleach_feed + rho_wet_ore

        vol_flowrate_initial = mass_flowrate_initial/rho_preleach

        if vol_flowrate_initial == 0:
            vol_flowrate_initial = wet_ore_mass_flowrate/rho_wet_ore
            tau = self.preleach_vol/vol_flowrate_initial
        else:
            tau = self.preleach_vol/vol_flowrate_initial

        mass_flowrate = mass_flowrate_inflow + \
                        math.exp(-time/tau) * (mass_flowrate_initial - mass_flowrate_inflow)

        tmp = self.preleach_phase.get_row(time)

        # Advance time
        time += self.time_step

        self.preleach_phase.add_row(time, tmp)

        self.preleach_phase.set_value('mass_flowrate', mass_flowrate, time)
        self.preleach_phase.set_value('mass_density', rho_preleach, time)

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
