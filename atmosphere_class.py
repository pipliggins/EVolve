"""
Definition of the Atmosphere class.

Contains data and methods relating to the running of FastChem and calculating
the composition of the atmosphere at each timestep.

"""

import os
import shutil
import pandas as pd
import numpy as np

from fastchem.python import pyfastchem
from fastchem_tools import saveChemistryOutput

import constants as cnst
import conversions as cnvs
import messages as msgs

class AtmosphereDef:
    """
    Stores the properties of the atmosphere + atmospheric results of each time step.
    
    Attributes:
        fastchem (bool): If True, use run EVo results through FastChem.
        fc_cheminput (dict): Stores the input for FC after conversion from EVo output.
        atmos_input (dict): Input to the atmosphere prior to mixing as a mol fraction from EVo.
        t_surf (float): Surface temperature (K).
        p_surf (float): Surface pressure (bar).
        mmw (float): Mean molecular mass of the atmosphere (g).
        composition (dict): Final composition for current and previous timestep (mol fraction).
        tracked_species (lst): The species which will be checked for atm. evolution.
        total_atomic_mass (dict): The total mass in the atmosphere for each element, updated after each eruption.
        height (float)     : atmospheric thickness (Km)

    """

    def __init__(self):
        self.fastchem = False
        self.fc_cheminput = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
        self.atmos_input = {}       
        self.t_surf = 288
        self.p_surf = 1
        self.mmw = None 
        self.composition = {}
        self.tracked_species = []
        self.total_atomic_mass = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
        self.height = 0     

    def atm_mass(self, planet):
        """
        Returns the mass of the atmosphere in g
        """
        return ((self.p_surf*1e8)*planet.surface_area)/planet.g    # p is converted from bar to Pa (kg m-1 s-2, *1e5) to g m-1 s-2 (*1e3)

    def init_store(self, planet):
        """
        Creates atmosphere_out.csv and saves initial atmospheric properties

        Creates a csv file to store the atmospheric composition and temperature,
        and instantiates with the starting values for the atmosphere.

        Args:
            planet (class): Active instance of the PlanetDef class
        
        Returns:
            None
        """

        df_atm = pd.DataFrame(self.composition, index=[-1])
        
        df_atm.insert(loc=0, column='t_step', value = 0)
        df_atm.insert(loc=1, column='age_yr', value = 0)
        df_atm.insert(loc=2, column='P_surf', value = self.p_surf)
        df_atm.insert(loc=3, column='height', value=self.scale_height(planet))
        df_atm.insert(loc=4, column='mmw', value = self.mmw)
        
        df_atm.to_csv('outputs/atmosphere_out.csv', sep='\t', index=False) 
        
        # delete the dataframe to reduce memory load
        del df_atm


    def initial_comp(self, planet, mantle, fastchem):
        """
        Uses the tracked species list to create the initial atmospheric composition.

        The tracked species list is compared to the provided atmospheric composition.
        Any species listed as tracked but not in the initial atmosphere is added as 0,
        and any species present in the initial atmospheric composition but not in the
        tracked species list is queried, with the option provided to add the species 
        to the tracked list or ignore and carry on the run.

        Args:
            planet (class): Active instance of the PlanetDef class
            mantle (obj): MantleDef object
            fastchem (bool): True if fastchem has been selected at the commandline
        
        Returns:
            fastchem (bool): enables selecting of FastChem to be used.
        """        

        allowable_names = self.tracked_species

        for name in allowable_names:
            if name in self.composition:
                self.composition[name] = [self.composition[name]]
            else:
                self.composition[name] = [0]

        names = [key for key in self.composition.keys() if key not in allowable_names]
        
        if names:
        
            answer = msgs.query_yes_no(f'The species {names} are listed in the initial atmospheric composition but are not being tracked to determine if atmospheric steady state has been reached.\nDo you wish to add these species to the tracked list?')

            if answer == True:
                for name in names:
                    self.tracked_species.append(name)
                    self.composition[name] = [0]
                print(f'Species {names} have been added to the tracked list.')
            
            else:
                answer = msgs.query_yes_no(f'Species {names} have NOT been added to the tracked list.\nSpecies being tracked are: {self.tracked_species}. Do you wish to continue?')

                if answer == False:
                    exit('Exiting EVolve.')

        if fastchem != True and (self.t_surf != mantle.t):
            print(f'Without FastChem enabled, the atmospheric composition will be for a gas at {mantle.t} K rather than the surface temperature {self.t_surf}.\n')
            
            answer = msgs.query_yes_no('Would you like to continue by turning on FastChem?', default='yes')

            if answer == True:
                print('Turning on FastChem option.')
                fastchem = True
            
            else:                        
                answer = msgs.query_yes_no(f'Do you wish to continue by setting the surface temperature to {mantle.t} K?', default=None)
                if answer == True:
                    print(f'Setting surface temperature to {mantle.t} K.')
                    self.t_surf = mantle.t
                else:
                    msgs.exiting('EVolve cannot generate an atmosphere cooler than magmatic temperature without FastChem.')
        
        self.init_store(planet)

        return fastchem

    
    def unitconvert_atomicmass2fc(self, masses=None):
        """
        Converts the atmospheric composition as atomic masses (weight of H, C, etc)
        to fastchem inputs, after H2 loss has occured.

        Args:
            masses (dict): optional - set of atomic MASSES (not fractions) to be input to FastChem
        
        Returns:
            fc (dict): The converted volatiles in FastChem units
        """

        if masses == None:
            mols = self.total_atomic_mass  # the total mass of each element in the atmosphere
        else:
            mols = masses
        
        # atm_norm = cnvs.mass2wt_frac(mols)  # mass fractions of each element in the atmosphere
        atm_norm = cnvs.wt2mol(cnvs.mass2wt_frac(mols))  # mass fractions of each element in the atmosphere
        
        # convert to fastchem units
        fc = {}
        for x, n in atm_norm.items():
            fc[x.lower()] = np.log10(n/atm_norm['H']) + 12.0
        
        return fc
   

    def fc_setup(self):
        """
        Edits the fastchem config file to that from the inputs folder.
        """

        files = ['parameters.dat', 'config.input']
        
        for f in files:
            if os.path.isfile(f'fastchem/input/{f}'):
                os.remove(f'fastchem/input/{f}')
        
            shutil.copyfile(f'inputs/{f}', f'fastchem/input/{f}')


    def fc_filecreate(self, t_step):
        """
        Creates the atmosphere and chemistry input files to run FastChem with.

        Creates atmosphere and chemcistry input files based on the chemistry input
        from EVo and the current temperature and pressure of the atmosphere.
        Stores these input parameters in a file in the outputs folder, fc_input.csv.

        Args:
            t_step (float): Current timestep

        Returns:
            None
        """
                
        fabun = open('fastchem/input/evolve_abundances.dat','w')

        fabun.write('#The chemistry input abundances from evo')
        fabun.write('\n')
        for ele in self.fc_cheminput.keys():
            fabun.write(f'{ele.upper()}  {self.fc_cheminput[ele]}')
            fabun.write('\n')
                    
        fabun.close()

        file_exists = os.path.isfile('outputs/fc_input.csv')
        
        if file_exists:
            with open('outputs/fc_input.csv', 'a') as f_inp:
                f_inp.write(f"{t_step} \t {self.t_surf:.3f} \t {self.p_surf:.16f} \t {self.fc_cheminput['h']} \t {self.fc_cheminput['c']} \t {self.fc_cheminput['o']} \t {self.fc_cheminput['s']} \t {self.fc_cheminput['n']}\n")
        
        else:
            with open('outputs/fc_input.csv', 'w') as f_inp:
                f_inp.write(f't_step \t T_K \t P_bar \t H \t C \t O \t S \t N\n')
                f_inp.write(f"{t_step} \t {self.t_surf:.3f} \t {self.p_surf:.16f} \t {self.fc_cheminput['h']} \t {self.fc_cheminput['c']} \t {self.fc_cheminput['o']} \t {self.fc_cheminput['s']} \t {self.fc_cheminput['n']}\n")
        
    
    def run_fastchem(self, planet, t_step):
        """
        Sets up FastChem input files, then runs FastChem.

        Args:
            planet (obj): PlanetDef object
            t_step (float): Current timestep
        
        Returns:
            None
        """

        # Create the abundance file
        self.fc_filecreate(t_step)

        # pyfastchem -----------------------
        fastchem = pyfastchem.FastChem('fastchem/input/parameters.dat', 1)

        input_data = pyfastchem.FastChemInput()
        output_data = pyfastchem.FastChemOutput()

        input_data.temperature = [float(self.t_surf)]
        input_data.pressure = [float(self.p_surf)]

        result = fastchem.calcDensities(input_data, output_data)

        # if convergence ok result == 0
        if result != 0:
            msgs.exiting(pyfastchem.FASTCHEM_MSG[result])       # This negates the need to check the terminal output for error codes.
        
        # -----------------------------------------

        self.fc_saveoutput(planet, t_step, output_data, fastchem)


    def fc_saveoutput(self, planet, tstep, output_data, fastchem):
        """
        Appends the whole FastChem output to fastchem_out.csv
        and saves results to atm.composition.
        
        Args:
            planet (class): Active instance of Planet class
            tstep (int): The current iteration timestep
            output_data (class): output from running fastchem
            fastchem (class): object from fastchem python run

        Returns:
            None
        """

        filename = 'fc_out.csv'
        
        mmw, tracked_sp = saveChemistryOutput(filename, self.t_surf, self.p_surf, output_data, fastchem, tstep, tracked_species=self.tracked_species)

        for sp in tracked_sp:
            try:
                self.composition[sp].append(tracked_sp[sp])
            except: # if tracked element isn't computed in evo/fastchem, e.g. He
                continue
        
        self.mmw = mmw
        self.height = self.scale_height(planet)

    def saveatmosphere(self, planet, new_comp, tstep):
        """
        Saves most recent iteration out to atmosphere_out.csv,
        updates the composition dict

        Args:
            planet (class): active instance of PlanetDef class
            new_comp (dict): The new, mixed atmospheric composition
            tstep (int): Timestep
        
        Returns:
            None
        """

        for sp in new_comp.keys():
            self.composition[sp].append(new_comp[sp])
            
        self.mmw = cnvs.calc_mmw(self.composition, i=-1)
        self.height = self.scale_height(planet)
        
        file_exists = os.path.isfile('outputs/atmosphere_out.csv')

        df_atm = pd.DataFrame(new_comp, index=[-1])
            
        df_atm.insert(loc=0, column='t_step', value = tstep)
        df_atm.insert(loc=1, column='age_yr', value = tstep*planet.timestep)
        df_atm.insert(loc=2, column='P_surf', value = self.p_surf)
        df_atm.insert(loc=3, column='height', value=self.height)
        df_atm.insert(loc=4, column='mmw', value = self.mmw)
    
        if file_exists:
            df_atm.to_csv('outputs/atmosphere_out.csv', mode='a', header=False, sep='\t', index=False)
        else:
            df_atm.to_csv('outputs/atmosphere_out.csv', sep='\t', index=False)
    
        # delete the dataframe to reduce memory load
        del df_atm

    def fc_saveatmosphere(self, planet, tstep):
        """
        Saves most recent iteration out to atmosphere_out.csv
        
        In the case that FastChem is used, the atmosphere_out.csv file needs saving to
        after FC is run, rather than within mix_atmosphere as happens when FC is ignored.

        Args:
            planet (class): active instance of PlanetDef class
            tstep (int): Timestep
        
        Returns:
            None
        """

        new_comp = {}

        for sp in self.composition.keys():
            new_comp[sp] = self.composition[sp][-1]
        
        file_exists = os.path.isfile('outputs/atmosphere_out.csv')

        df_atm = pd.DataFrame(new_comp, index=[-1])
                
        df_atm.insert(loc=0, column='t_step', value = tstep)
        df_atm.insert(loc=1, column='age_yr', value = tstep*planet.timestep)
        df_atm.insert(loc=2, column='P_surf', value = self.p_surf)
        df_atm.insert(loc=3, column='height', value=self.scale_height(planet))
        df_atm.insert(loc=4, column='mmw', value = self.mmw)
        
        if file_exists:
            df_atm.to_csv('outputs/atmosphere_out.csv', mode='a', header=False, sep='\t', index=False)
        else:
            df_atm.to_csv('outputs/atmosphere_out.csv', sep='\t', index=False) 
        
        # delete the dataframe to reduce memory load
        del df_atm
    
    def update_psurf(self, planet, volc):
        """
        Calculates the new surface pressure after a volcanic degassing step.

        Calculates the new surface pressure by using the mass of the outgassed component
        to calculate it's partial pressure and therefore the new surface pressure.

        Args:
            planet (class): The active instance of the Planet class
            volc (class): The active instance of the Volcano class

        Returns:
            p_new (float): The new surface pressure after volcanic outgassing (bar)
            p_volc (float): The surface pressure of the erupted volcanic gases only
        """

        A = planet.surface_area # Surface area of the silicate body, m^2

        m_volc = (volc.melt_mass * volc.WgT)             # this should be grams

        p_volc = ((m_volc * planet.g)/A)*1e-8            # convert g m-1 s-2 to bar

        p_new = self.p_surf + p_volc   
        
        return p_new, p_volc

    def update_store(self, volc, melt_mass):
        """
        Updates the store of total atomic masses in the atmosphere with the
        contents of the volcanic gas.

        args:
            volc (obj): VolcanoDef object
            melt_mass (float): the mass of the melt + gas system erupted in EVo
        
        Returns:
            None
        """

        # append the erupted gas atomic masses to atmosphere store.
        for key in self.total_atomic_mass.keys():
            wt = volc.atomicM_gas[key.lower()]*1e-6*melt_mass       # atomicM_gas is still the weight fractions of the system, just only the ones in the gas phase.
            self.total_atomic_mass[key] += wt # mass of each element in the atmosphere (circumventing fact we only track a few species from fastchem)

    
    def mix_atmosphere(self, p_new, p_volc):
        """
        Calculates the chemistry of the new atmosphere after volcanic degassing.

        Combined the pre-existing atmospheric composition with the volcanic input to
        generate the new atmospheric composition as a mol fraction.

        Args:
            p_new (float): The new atmospheric pressure after degassing
            p_volc (float): The surface pressure of the erupted volcanic gases only
        
        Returns:
            new_comp (dict): the new atmospheric composition after evo has erupted
            p_new (float): the new surface pressure
        """

        new_comp = {}  # stores the new atmospheric composition as mol fractions ready to save in csv file.
        
        for sp in self.composition:
            if sp in self.atmos_input:
                pi_volc = self.atmos_input[sp]*p_volc   # atmos_input is the input from evo/fc as mol fractions
                pi_init = self.composition[sp][-1]*self.p_surf
                pi_new = pi_volc + pi_init
                new_mol_frac = pi_new/p_new
                new_comp[sp] = new_mol_frac
            else:
                new_comp[sp] = 0.0
                
        return new_comp, p_new
    
    def scale_height(self, planet):
        """
        Calculates the atmospheric scale height (km)

        Following method of Ortenzi 2020, adjusted so that t_atm is
        the specified surface temperature.

        Args:
            planet (class): Active instance of the Planet class

        Returns:
            H (float): Atmoshperic scale height (km)
        """

        H = (cnst.R * self.t_surf)/((self.mmw / 1e3) * planet.g)        # pressure scale height in m, mmw converted to kg/mol

        return H/1e3        # scale height, in km
    
    
    def check_evolution(self, tol=1e-3):
        """
        Check if the atmospheric composition is still changing.

        Compares the new atmospheric composition to the previous step, to establish rate
        of change. Returns True if the atmospheric composition is still changing at a 
        resolution greater than tol.

        Args:
            tol (float): The tolerance on the rate of change
        
        Returns:
            bool: True if still changing, False if difference < tol.
        """
        
        for ele in self.composition.keys():
            if abs(self.composition[ele][0]-self.composition[ele][1])/self.composition[ele][0] > tol:
                # delete the first element of each species in composition so it just reflects the current composition
                for ele in self.composition.keys():
                    del self.composition[ele][0]
                return True

        # If none of the elements have changed by more than tol, trigguring true, then
        return False