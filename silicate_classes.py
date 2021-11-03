"""
Definitions of the Planet, Mantle and Volcano classes.

Stores class definitions for use in EVolve.
"""

import os, csv, ruamel.yaml
import subprocess as sp
from subprocess import PIPE
import numpy as np
import pandas as pd

import constants as cnst
import conversions as cnvs
import messages as msgs
from evo.density import den_calc_spera2000 as density

class PlanetDef:
    """
    Stores the properties of the planet and temporal evolution as a whole. 
    
    Attributes:
        age (int):             Age of the planet at time 0, or start of degassing (years)
        core_mass (float):     Mass of the rocky part of the planet; (core, mantle + crust). (grams).
        radius (float):        Radius of the rocky part of the planet
        g (float):             gravitational acceleration at the surface
        melt_flux (float):     Volume of melt erupted to the surface, m3/yr
        set_melt_atomic_weights (bool): If true, rather than setting H2O, CO2 etc in evo,
                                they are converted to H, C, S, N and input to evo to 
                                calculate true speciation based on the fO2.
        irradiation (float):   Stellar irradiance
        timestep (int):        Timestep for each model iteration, years
        max_timestep (int):    Max number of timesteps to execute
        snapshot (bool):       If true, only do one timestep.
        h_loss_rate (float):    If known, the rate of H loss from the atmosphere in g s-1.
        total_volatiles (dict):The total MASS of each element in the PLANET. 
        surface_area:          Returns the surface area of the planet.
    """

    def __init__(self):
        self.age = 0
        self.core_mass = 0
        self.radius = 0
        self.g = 9.81
        self.melt_flux = 30.0e9
        self.set_melt_atomic_weights = 'true'
        self.irradiation = 0
        self.timestep = 1000
        self.max_timestep = 200
        self.snapshot = False
        self.h_loss_rate = None
        self.total_volatiles = {}

    @property
    def surface_area(self):
        return 4 * np.pi * self.radius **2

    def setup_conservation(self, mantle, atm):
        """
        Calculates the total mass of volatiles in the whole planetary system, to allow
        checks later for volatile loss somewhere.

        Args:
            mantle (class): active instance of MantleDef class
            atm (class): active instance of Atmosphere class

        Returns:
            None
        """

        mantle_h, mantle_c, mantle_s, mantle_n = mantle.calc_atomic_masses(mantle.mantle_volatiles)  # atomic fractions in ppm

        mantle_h = mantle_h*1e-6*mantle.mass
        mantle_c = mantle_c*1e-6*mantle.mass
        mantle_s = mantle_s*1e-6*mantle.mass
        mantle_n = mantle_n*1e-6*mantle.mass

        atm_mass = atm.atm_mass(self)   # g

        atm_wt_frac = cnvs.mol2wt(atm.composition, i=0)
        atm_comp_mass = {v:atm_wt_frac[v]*atm_mass for v in atm.composition.keys()} # mass of each molecule in atmosphere (weight fraction * atmospheric mass)

        H = 2*cnst.m['h']*((atm_comp_mass['H2O']/cnst.m['h2o']) + (atm_comp_mass['H2']/cnst.m['h2']) + (atm_comp_mass['H2S']/cnst.m['h2s']) + 2*(atm_comp_mass['CH4']/cnst.m['ch4']) + (atm_comp_mass['C2H2']/cnst.m['c2h2'])) + cnst.m['h']*((atm_comp_mass['HCN']/cnst.m['hcn']) + 3*(atm_comp_mass['NH3']/cnst.m['nh3']))

        C = cnst.m['c'] *((atm_comp_mass['CO2']/cnst.m['co2']) + (atm_comp_mass['CO']/cnst.m['co']) + (atm_comp_mass['CH4']/cnst.m['ch4']) + (atm_comp_mass['HCN']/cnst.m['hcn']) + 2*(atm_comp_mass['C2H2']/cnst.m['c2h2']))

        S = cnst.m['s']*((atm_comp_mass['SO2']/cnst.m['so2']) + 2*(atm_comp_mass['S2']/cnst.m['s2']) + (atm_comp_mass['H2S']/cnst.m['h2s']))

        N = cnst.m['n'] * (2*(atm_comp_mass['N2']/cnst.m['n2']) +  (atm_comp_mass['HCN']/cnst.m['hcn']) + (atm_comp_mass['NH3']/cnst.m['nh3']))

        O = cnst.m['o']*((atm_comp_mass['H2O']/cnst.m['h2o']) + 2*(atm_comp_mass['O2']/cnst.m['o2']) + 2*(atm_comp_mass['CO2']/cnst.m['co2']) + (atm_comp_mass['CO']/cnst.m['co']) + 2*(atm_comp_mass['SO2']/cnst.m['so2']))

        # mass of each element in the atmosphere
        atm.total_atomic_mass['H'] += H
        atm.total_atomic_mass['C'] += C
        atm.total_atomic_mass['S'] += S
        atm.total_atomic_mass['N'] += N
        atm.total_atomic_mass['O'] += O
        
        self.total_volatiles['h'] = H + mantle_h
        self.total_volatiles['c'] = C + mantle_c
        self.total_volatiles['s'] = S + mantle_s
        self.total_volatiles['n'] = N + mantle_n
    
    def check_conservation(self, mantle, atm, nocrust=False):
        """
        Checks that volatiles are being conserved within the system
        (except H if hydrogen escape is on).

        Args:
            mantle (class): Active instance of MantleDef class
            atm (class): Active instance of PlanetDef class
            nocrust (bool): If True, crust formation does not occur

        Returns:
            None
        """

        mantle_vols = {}
        crust_vols_mass = {}

        vols = ['H2O', 'CO2', 'N', 'S']

        for v in vols:
            mantle_vols[v] =  mantle.mantle_volatiles[v]/1e6
            mantle_vols[v] = (mantle_vols[v]*mantle.mass)/cnst.m[v.lower()] # number of mols of each species (mass/molar_mass)

            if nocrust == False:
                crust_vols_mass[v] = mantle.crust[v]/1e6
                crust_vols_mass[v] = (crust_vols_mass[v]*mantle.crust['mass'])/cnst.m[v.lower()]
            else:
                crust_vols_mass[v] = 0

        H = 2*cnst.m['h']*(mantle_vols['H2O'] + crust_vols_mass['H2O']) + atm.total_atomic_mass['H']

        C = cnst.m['c'] * (mantle_vols['CO2'] + crust_vols_mass['CO2']) + atm.total_atomic_mass['C']

        S = cnst.m['s'] * (mantle_vols['S'] + crust_vols_mass['S']) + atm.total_atomic_mass['S']

        N = cnst.m['n'] * (mantle_vols['N'] + crust_vols_mass['N']) + atm.total_atomic_mass['N']

        if abs(self.total_volatiles['h'] - H)/H > 1e-10 or abs(self.total_volatiles['c'] - C)/C > 1e-10 or abs(self.total_volatiles['s'] - S)/S > 1e-10 or abs(self.total_volatiles['n'] - N)/N > 1e-10:
            exit('Error: volatiles are not being adequately conserved.')


class MantleDef:
    """
    Stores the properties of the mantle. Holds the results of each time step.
    
    Holds methods for calculating the size and volatile content of a melt packet,
    set up the input files for EVo in accordance with the properties of this melt
    packet, and write the properties (vol content, mass etc) of the mantle after
    melting at each timestep out to a csv file.

    Attributes:
        t (float): Mantle temperature (K); starting T for EVo
        mass (float): Mass of the mantle
        fo2_buffer (str): The rock buffer fo2 will be quoted as relative to
        fo2 (float): Mantle fO2 relative to the rock buffer
        ie_ratio (float): Intrusive:extrusive ratio for magma produced in the mantle
        F (float): Mantle averaged melt fraction for batch melting
        composition (dict): dry major oxide composition melt produced by mantle melting
        mantle_volatiles (dict): Volatile content of the mantle after each timestep (ppm)
        crust (dict): volatile concs that didn't degass according to IE ratio + crustal mass

    """

    def __init__(self):
        self.t = 1473.25
        self.mass = 1e10
        self.fo2_buffer = 'IW'
        self.fo2 = 0
        self.ie_ratio = 0.1  # Fraction of melt which reaches the surface to degass. Use as a volcanic flux? - see spaargaren, tosi & others
        self.F = 0.1
        self.composition = None
        self.mantle_volatiles = {}
        self.crust = {}                 # store of volatiles that didn't degass according to IE ratio PL: return these to mantle or store as inaccessible?
    
    def init_store(self):
        """
        Creates mantle_out.csv and saves initial mantle properties

        Creates a csv file to store the mantle and crustal masses + volatile contents,
        and instantiates with the starting values for the mantle.

        Args:
            None
        
        Returns:
            None
        """

        data = {'t_step':0, 'age_yr':0, 'mantle_mass':'{:3.4e}'.format(self.mass), self.fo2_buffer:self.fo2}

        for key in self.mantle_volatiles.keys():
            data['mantle_'+key] = '{:.4f}'.format(self.mantle_volatiles[key])

        data['melt_H2O'], data['melt_CO2'], data['melt_S'], data['melt_N'] = 0,0,0,0
        
        data['crust_mass'] = 0

        data['crust_thickness'] = 0

        for key in self.mantle_volatiles.keys():
            data['crust_'+key] = 0
        
        with open('outputs/mantle_out.csv', 'w') as f:
            w = csv.DictWriter(f, data.keys(), delimiter='\t')
            w.writeheader()
            w.writerow(data)


    def batch_melt(self, F, vol):
        """
        Finds the concentration of a volatile in the melt phase.

        Executes the batch melting equation for each volatile in the system to calculate
        the concentration of each volatile in the melt phase, in ppm.

        Args:
            F (float): the melt fraction of the system
            vol (str): the volatile name
        
        Returns:
            Concentration of volatile 'vol' in the melt phase in ppm
        """

        if vol.lower() != 'n':
            return self.mantle_volatiles[vol]/(cnst.D[vol.lower()] + F*(1 - cnst.D[vol.lower()]))
        else:       #N
            return self.mantle_volatiles[vol]/(cnst.D['n_fefeo'] + F*(1 - cnst.D['n_fefeo']))  # PL: Only using the reduced form at the moment, either turn into linear equation or split according to which is closer.


    def mix_reserviors(self, m1, x1, m2, x2, removing=True):
        """
        Calculate the new volatile concentration of a reservoir (Tosi et al 2017, eq22).

        The volatile concentration of a melt/solid reservoir is calculated according to
        the mass and concentration of a volatile in each of the 2 original reservoirs,
        depending on whether one is being added to or removed from the other.

        Args:
            m1 (float): mass of the initial reservior..
            x1 (float): concentration (ppm) of the volatile in the initial reservior
            m2 (float): mass of the body being added/removed
            x2 (float): concentration (ppm) in the body being added/removed
            removing (bool): If True, mass will be removed. If false, it will be added.
        
        Returns:
            x_new (float): concentration (ppm) in the new reservoir
        """

        if removing == True:
            x_new = (x1*m1 - x2*m2)/(m1-m2)
        
        else:
            x_new = (x1*m1 + x2*m2)/(m1+m2)
        
        return x_new

    def generate_melt(self, planet, atm, tstep, nocrust=False):
        """ 
        Generates the mass and volatile content of a melt that will be 'erupted' in EVo.

        Generates a melt phase with size based on the amount of mantle melting and
        the intrusive to extrusive ratio. Volatile contents are controlled by the 
        batch melting equation and the partition coefficients of each species.

        Updates the mantle volatile storage to remove the volatiles which are
        now present in the melt phase and are being erupted/ stored in the crust.

        Stores the volatiles and melt mass remaining as the intrusive component into the
        'crust' dictionary - these remain in the crust and are not returned to the mantle
        or erupted to the surface.
        
        Args:
            planet (class): The active instance of the planet class
            atm (class): The active instance of the Atmosphere class
            tstep (int): The current iteration timestep
            nocrust (bool): If True, no crust reservoir is formed.
        
        Returns:
            melt_volatile (dict): the volatile concentrations in the melt, in ppm
            melt_mass (float): The mass of the melt body being erupted to the surface by EVo
        """
        
        # Calculate mass of melt (silicate only) based on volume + density at surface pressure
        melt_density = density(cnvs.wt2mol(self.composition), self.t, atm.p_surf)*1e3  # g/m3
        melt_mass = melt_density*planet.melt_flux # mass of melt erupted to the surface in one year

        melt_mass *= planet.timestep  # the amount of melt erupted to the surface (volcanic flux), g, in the whole time period.
        total_melt_mass = melt_mass * 1/self.ie_ratio   # total mass of melt produced in the mantle
        
        intrusive_mass = total_melt_mass - melt_mass  # amount of that melt which will get stored in the crust undegassed

        melt_volatile = {}

        
        for species in self.mantle_volatiles.keys():
            melt_volatile[species] = self.batch_melt(self.F, species) # CONCENTRATION (ppm) of volatile in the melt.

            if nocrust == True:
                # Update mantle concs by removing ONLY species in the melt, so implicitly the residual mush stays in the mantle.
                # In this case only remove as much melt as is going into the volcanic system. all other would go into crust -> return to mantle, so just cutting out middle man.
                self.mantle_volatiles[species] = self.mix_reserviors(self.mass, self.mantle_volatiles[species], melt_mass, melt_volatile[species], removing=True)

            elif nocrust == False:            
                # Update mantle concs by removing ONLY species in the melt, so implicitly the residual mush stays in the mantle.
                # in this case both intrusive & extrusive melt removed from mantle, the intrusive will go into the crust, extrusive to volcano.
                self.mantle_volatiles[species] = self.mix_reserviors(self.mass, self.mantle_volatiles[species], total_melt_mass, melt_volatile[species], removing=True)
                
                # Update crustal store
                if species in self.crust:
                    self.crust[species] = self.mix_reserviors(self.crust['mass'], self.crust[species], intrusive_mass, melt_volatile[species], removing=False)
                else:
                    self.crust[species] = melt_volatile[species]  # when there is no pre-existing crust (in the model), the CONCENTRATION of volatile in crust and melt will be equal.
                
        if nocrust == True:
            self.mass -= melt_mass      # only remove the stuff that's going into the volcano.
        
        elif nocrust == False:
            # Update the mantle mass
            self.mass -= total_melt_mass
        
            if 'mass' in self.crust:
                self.crust['mass'] += intrusive_mass
            else:
                self.crust['mass'] = intrusive_mass
        
        return melt_volatile, melt_mass
    
    def update_crust(self, planet, volc, atm, concs, WgT, melt_mass, melt_v, tstep, nocrust=False):
        """
        Updates the crustal store with the degassed melt left after an EVo run.

        Adds the degassed melt back into the crust or mantle, then saves out the mantle results, plus saves the total
        atomic volatiles to the atmosphere class.

        Args:
            planet (class): The active instance of the Planet class
            volc (class): The active instance of the VolcanoDef class
            concs (dict): The concentrations of volatiles left over in the melt as wt fracs.
            WgT (float): gas weight fraction
            melt_mass (float): Undegassed mass of the melt before it went through EVo.
            melt_v (dict): The concentrations of volatiles in the melt setting up EVo
            tstep (int): The current iteration timestep
            nocrust (bool): If True, degassed melt added straight back into mantle

        Returns:
            None (updates the mantle.crust dictionary internally)
        """

        degassed_mass = melt_mass*(1-WgT)

        h, c, s, n = volc.atomicM_melt['h'], volc.atomicM_melt['c'], volc.atomicM_melt['s'], volc.atomicM_melt['n']  # volatiles in the degassed melt converted into atomic weight fractions (ppm)

        def system2melt(v, melt_mass, degassed_mass):
            v = v*1e-6  # ppm to wt frac
            v = v * melt_mass
            v = v/degassed_mass
            v = v*1e6  # wt frac to ppm
            return v
        
        h = system2melt(h, melt_mass, degassed_mass)    # ppm in the remaining melt
        c = system2melt(c, melt_mass, degassed_mass)
        s = system2melt(s, melt_mass, degassed_mass)
        n = system2melt(n, melt_mass, degassed_mass)

        h2o, co2, s, n = cnvs.atomic_mass2mantle_volatiles([h, c, s, n])

        converted_melt = {}
        converted_melt['H2O'] = h2o
        converted_melt['CO2'] = co2
        converted_melt['S'] = s
        converted_melt['N'] = n


        if nocrust == True:
            # add the degassed melt and it's volatiles back into the mantle.
            for species in converted_melt.keys():
                self.mantle_volatiles[species] = self.mix_reserviors(self.mass, self.mantle_volatiles[species], degassed_mass, converted_melt[species], removing=False)

                self.crust[species] = 0

            self.mass += degassed_mass
            self.crust['mass'] = 0
            self.crust['thickness'] = 0
        
        elif nocrust == False:
            # add degassed melt and it's volatiles into the crust
            for species in converted_melt.keys():
                self.crust[species] = self.mix_reserviors(self.crust['mass'], self.crust[species], degassed_mass, converted_melt[species], removing=False)
            
            self.crust['mass'] += degassed_mass
            
            crust_density = density(cnvs.wt2mol(self.composition), atm.t_surf, atm.p_surf)*1e3  # approx crustal density at surface PT conditions

            crust_vol = self.crust['mass']/crust_density

            thickness = crust_vol/planet.surface_area   # PL: this assumes a flat surface not a sphere...

            self.crust['thickness'] = thickness

        self.save_mantle(planet, melt_v, tstep)    

    def calc_atomic_masses(self, melt_vols):
        """
        Takes melt volatile contents as ppm H2O, CO2 etc and calculates the atomic
        weight fractions of H, S, C and N.

        Args:
            melt_vols (dict): Containing the melt volatile contents in ppm

        Returns:
            tuple: (H, C, S, N) as atomic mass fractions in ppm.
        """

        h2o = melt_vols['H2O']
        co2 = melt_vols['CO2']
        s = melt_vols['S'] # These are already == to atomic mass fractions as they dissolve as ions with no O.
        n = melt_vols['N']

        if 'H2' in melt_vols:
            h2 = melt_vols['H2']
        else:
            h2 = 0
        
        if 'CO' in melt_vols:
            co = melt_vols['CO']
        else:
            co = 0
        
        if 'CH4' in melt_vols:
            ch4 = melt_vols['CH4']
        else:
            ch4 = 0

        h = 2*cnst.m['h']*((h2o/cnst.m['h2o']) + (h2/cnst.m['h2']) + 2*(ch4/cnst.m['ch4']))
        c = cnst.m['c']*((co2/cnst.m['co2']) + (co/cnst.m['co']) + (ch4/cnst.m['ch4']))

        return h, c, s, n
    
    def set_evo_env(self, volc, atm, set_atomic_weights=False):
        """
        Edits the env.yaml input file for EVo.
        
        Updates env.yaml to match the current properties of the mantle derived melt.
        Sets the melt fO2 to be equal to mantle fO2, and the final surface pressure
        according to the atmospheric properties.

        Args:
            volc (class): The active instance of the VolcanoDef class
            atm (class): The active instance of the Atmosphere class
            set_atomic_weights (bool): If true, Evo is provided with atomic weight fractions to 
                                find the speciation & volatile saturation point independantly.
        
        Returns:
            None
        """

        file_loc = volc.env
        melt_vols = volc.melt_volatiles
        
        yaml = ruamel.yaml.YAML()   # ensures order of keys, comments and whitespace of file are maintained.

        with open(file_loc) as f:
            env = yaml.load(f)
        
        # calculate the system atomic mass fractions, for storage & use if atomic_set_weights == True.
        h, c, s, n = self.calc_atomic_masses(melt_vols)
        volc.atomicM_sys['h'] = h
        volc.atomicM_sys['c'] = c
        volc.atomicM_sys['s'] = s
        volc.atomicM_sys['n'] = n
        
        for prop in env:
            if prop == 'T_START':
                env[prop] = float(self.t)
                continue
            
            elif prop == 'P_STOP':
                env[prop] = cnvs.round_proper(float(atm.p_surf), 3)
                continue

            elif prop == 'FO2_buffer':
                env[prop] = self.fo2_buffer
                continue
            
            elif prop == 'FO2_buffer_START':
                env[prop] = float(self.fo2)
                continue

            elif prop == 'CO_MODEL':
                if self.fo2 < 1:
                    env[prop] = 'armstrong2015'
                else:
                    env[prop] = 'None'

            elif prop == 'CH4_MODEL':
                if self.fo2 < 1:
                    env[prop] = 'ardia2013'
                else:
                    env[prop] = 'None'
            
            if set_atomic_weights == False:
            
                if prop == 'FIND_SATURATION':
                    env[prop] = True
                    continue
                
                elif prop == 'ATOMIC_MASS_SET':
                    env[prop] = False
                    continue
                
                elif prop == 'WTH2O_SET':
                    env[prop] = True
                    continue
                
                elif prop == 'WTH2O_START':
                    env[prop] = float(melt_vols['H2O']/1e6)     # convert ppm to weight fraction
                    continue
            
                elif prop == 'WTCO2_SET':
                    env[prop] = True
                    continue
                
                elif prop == 'WTCO2_START':
                    env[prop] = float(melt_vols['CO2']/1e6)
                    continue

                elif prop == 'SULFUR_SET':
                    env[prop] = True
                    continue

                elif prop == 'SULFUR_START':
                    env[prop] = float(melt_vols['S']/1e6)
                    continue

                elif prop == 'NITROGEN_SET':
                    env[prop] = True
                    continue
                
                elif prop == 'NITROGEN_START':
                    env[prop] = float(melt_vols['N']/1e6)
                    break
            
            elif set_atomic_weights == True:

                if prop == 'FIND_SATURATION':
                    env[prop] = False
                    continue

                elif prop == 'ATOMIC_MASS_SET':
                    env[prop] = True
                    continue 
                
                elif prop == 'ATOMIC_H':
                    env[prop] = float(h)
                    continue

                elif prop == 'ATOMIC_C':
                    env[prop] = float(c)
                    continue

                elif prop == 'ATOMIC_S':
                    env[prop] = float(s)
                    continue

                elif prop == 'ATOMIC_N':
                    env[prop] = float(n)
                    continue
                
                elif prop == 'WTH2O_SET':
                    env[prop] = False
                    continue
                 
                elif prop == 'WTCO2_SET':
                    env[prop] = False
                    continue

                elif prop == 'SULFUR_SET':
                    env[prop] = False
                    continue
                
                elif prop == 'NITROGEN_SET':
                    env[prop] = False
                    break

            
        with open(file_loc, 'w') as f:
            yaml.dump(env, f)

    def save_mantle(self, planet, melt_v, tstep):
        """
        Save mantle and crustal properties out to a csv file.
        
        Saves the current mantle mass, fO2 and volatile content, along with crustal mass
        and volatile content, out to a csv file for each timestep.

        Args:
            planet (class): The active instance of the planet class
            melt_v (dict): The volatile content of the melt about to be input to EVo
            tstep (int): The current iteration timestep
        Returns:
            None
        """

        data = {'t_step':tstep, 'age_yr':tstep*planet.timestep, 'mantle_mass':'{:3.4e}'.format(self.mass), self.fo2_buffer:self.fo2}

        for key in self.mantle_volatiles.keys():
            data['mantle_'+key] = '{:.4f}'.format(self.mantle_volatiles[key])

        for key in melt_v.keys():
            data['melt_'+key] = '{:.4f}'.format(melt_v[key])
        
        data['crust_mass'] = '{:3.4e}'.format(self.crust['mass'])
        data['crust_thickness'] = '{:3.4e}'.format(self.crust['thickness'])

        for key in self.mantle_volatiles.keys():
            data['crust_'+key] = '{:.4f}'.format(self.crust[key])        
        
        with open('outputs/mantle_out.csv', mode='a') as f:
            w = csv.DictWriter(f, data.keys(), delimiter='\t')
            w.writerow(data)
    
class VolcanoDef:
    """
    Stores properties + methods for the volcanic system, relating to the EVo model.

    Holds methods for running and pulling data from the output of EVo.

    Attributes:
        atm (class): Active instance of the Atmosphere class.
        
        env (str):  Locations of input files for EVo as dictated by the input.yaml file.
        chem (str): ""
        out (str):  ""

        melt_mass (float):     Mass of the melt packet being erupted by EVo (g).
        melt_volatiles (dict): Volatile content of the melt at depth.

        atomicM_sys (dict): Atomic masses of whole system (ppm).
        atomicM_gas (dict): Atomic mass fractions of the gas phase at the final step (ppm).
        atomicM_melt (dict): Atomic mass fractions of the melt phase at the final step (ppm).
        mol_fracs (dict):   Speciation of the final gas phase as mol fraction.
        wt_fracs (dict):    Speciation of the final gas phase as weight fraction.
        melt_fracs (dict):  Volatiles left in the final melt as weight fraction.
        mmw (float):        Mean molecular weight of the final gas phase.
        WgT (float):        Final weight fraction of the gas phase.
        init_graph_sat (bool): Records whether melt was graphite saturated at volatile saturation.
    
    """

    def __init__(self, atm):
        self.atm = atm 
        
        self.env = "environment_file_location"
        self.chem = "chemistry_file_location"
        self.out = "output_file_location"
        
        self.melt_mass = 0
        self.melt_volatiles = {}
        
        self.atomicM_sys = {'c':0, 'h':0, 'n':0, 'o':0, 's':0}
        self.atomicM_gas = {'c':0, 'h':0, 'n':0, 'o':0, 's':0}      # SYSTEM MASS FRACTIONS (PPM) of volatiles only in the GAS phase.
        self.atomicM_melt = {'c':0, 'h':0, 'n':0, 'o':0, 's':0}      # SYSTEM MASS FRACTIONS (PPM) of volatiles only in the MELT phase.
        self.mol_fracs = {}
        self.wt_fracs = {}
        self.melt_fracs = {}
        self.mmw = 0
        self.WgT = 0
        self.init_graph_sat = False


    
    def get_atomicM_gas(self, WgT):
        """
        Calculates the atomic masses of elements in the gas phase only from EVo's output.

        Calculates the gas phase atomic masses in ppm from the EVo dataframe, then returns
        in dictionary format.

        Args:
            WgT (float): The gas mass fraction from the EVo run.

        Returns:
            atomicM_gas (dict): A dictionary of all the gas phase's volatile atomic mass fractions
        """

        mf = self.mol_fracs
        mjMj = [mf[ele] * cnst.m[ele.lower()] for ele in mf]

        self.atomicM_gas['o'] = cnst.m['o'] * 1e6 * ((WgT / sum(mjMj)) * (mf['H2O'] + 2*mf['O2'] + 2*mf['CO2'] + mf['CO'] + 2*mf['SO2']))

        self.atomicM_gas['h'] = 2 * cnst.m['h'] * 1e6 * ((WgT / sum(mjMj)) * (mf['H2O'] + mf['H2'] + 2*mf['CH4'] + mf['H2S']))

        self.atomicM_gas['c'] = cnst.m['c'] * 1e6 * ((WgT / sum(mjMj)) * (mf['CO'] + mf['CO2'] + mf['CH4']))

        self.atomicM_gas['s'] = cnst.m['s'] * 1e6 * ((WgT / sum(mjMj)) * (2*mf['S2'] + mf['H2S'] + mf['SO2']))

        self.atomicM_gas['n'] = cnst.m['n'] * 1e6 * ((WgT / sum(mjMj)) * (2*mf['N2']))

        return self.atomicM_gas

    
    def get_molmelt_fracs(self, df):
        """
        Retrieves the gas phase composition and melt volatile contents from EVo output.

        Finds the final step gas phase composition as both mol and weight fractions,
        and the melt volatile contents as weight fractions, and returns as dicts.

        Args:
            df (dataframe): The EVo output file opened as a Pandas dataframe

        Returns:
            mol_fracs (dict): final EVo gas composition as mol fractions
            wt_fracs (dict): final EVo gas composition as weight fractions
            melt_fracs (dict): final EVo melt volatile contents as weight fractions
        """

        lst_mol = ['mH2O', 'mH2', 'mO2', 'mCO2', 'mCO', 'mCH4', 'mSO2', 'mS2', 'mH2S', 'mN2']
        lst_wt = ['wH2O', 'wH2', 'wO2', 'wCO2', 'wCO', 'wCH4', 'wSO2', 'wS2', 'wH2S', 'wN2']
        lst_melt = ['H2O_melt', 'H2_melt', 'CO2_melt', 'CO_melt', 'CH4_melt', 'Stot_melt','N_melt']

        for elem, elew in zip(lst_mol, lst_wt):
            name = elem.lstrip('m')
            self.mol_fracs[name] = df[elem].iloc[-1]
            self.wt_fracs[name] = df[elew].iloc[-1]

        for ele in lst_melt:
            if ele != 'Stot_melt':
                name = ele.rstrip('_melt')  
            else:
                name = 'S'
            self.melt_fracs[name] = df[ele].iloc[-1]
        
        return self.mol_fracs, self.wt_fracs, self.melt_fracs
    
    def run_evo(self, tstep):
        """
        Runs EVo and saves out the final pressure step to volc_out.csv

        Runs the volcanic degassing model EVo in a shell, then stores the data in the Volcano class,
        and appends the final line of the EVo output file in the volc_out.csv file.

        Args:
            tstep (int): The current iteration timestep

        Returns:
            None
        """

        # Make the input file locations relative to the EVo file
        env = "../"+self.env
        chem = "../"+self.chem

        if self.out != 'None':
            out = "../"+self.out
            try:
                sp.run(f"cd evo && python dgs.py {chem} {env} --output {out}", shell=True, check=True, stderr=PIPE) # Change directory to EVo and run model.
            except sp.CalledProcessError as error:
                out =  error.stderr.decode("utf-8")
                out_split = out.split('\n')
                msgs.exiting(out_split[-2] + out_split[-1] + f'\nExiting: EVo failed to run at timestep {tstep}. :(')
            
        else:
            try:
                sp.run(f"cd evo && python dgs.py {chem} {env}", shell=True, check=True)
            except Exception as e:
                msgs.exiting(f'Exiting: EVo failed to run. :(\nError reported: {e}')
        
        with open("evo/Output/dgs_output.csv", 'r') as evo_output:
            evo_df = pd.read_csv(evo_output, sep='\t', comment='#')
        
        # Check the planet surface pressure and final EVo run pressure match

        # PL: UNCOMMENT THIS ONCE TESTING IS DONE!!!!!!!!!!!!!!!!!!!!!
        # if evo_df['P'].iloc[-1] != self.atm.p_surf:
        #     msgs.exiting('Error: Final pressure of EVo run does not match required surface pressure.')
        
        self.get_molmelt_fracs(evo_df)
        self.mmw = evo_df['mol_mass'].iloc[-1]
        self.WgT = evo_df['Gas_wt'].iloc[-1] / 100  # convert to weight fraction from percentage
        self.init_graph_sat = True if evo_df['graph_melt'].iloc[0] != 0.0 else False
        self.get_atomicM_gas(self.WgT)      # save the atomic weight fractions (ppm) in the gas phase only.

        for k in self.atomicM_gas.keys():
            self.atomicM_melt[k] = self.atomicM_sys[k] - self.atomicM_gas[k]    # rest of the system volatiles must be in the melt

        # Save bottom line to a csv file
        file_exists = os.path.isfile('outputs/volc_out.csv')

        if file_exists:
            evo_df = evo_df.iloc[[-1]]
            evo_df.insert(loc=0, column='t_step', value = tstep)
            evo_df.insert(loc=1, column='init_graph_sat', value = self.init_graph_sat)
            evo_df.to_csv('outputs/volc_out.csv', mode='a', header=False, sep='\t', index=False)
        else:
            evo_df = evo_df.iloc[[-1]]
            evo_df.insert(loc=0, column='t_step', value = tstep)
            evo_df.insert(loc=1, column='init_graph_sat', value = self.init_graph_sat)
            evo_df.to_csv('outputs/volc_out.csv', sep='\t', index=False) 
        
        # delete the dataframe to reduce memory load - necessary info now stored in volc + the csv
        del evo_df

