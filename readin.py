"""
Sets up EVolve's run parameters using the relevent YAML input files.

Creates and sets the parameters of each class by reading in the relevant
YAML files. Also stores all the initial conditions in a text file in the
outputs folder.
"""

import yaml
import inspect
import sys
import os
import ruamel.yaml

ryaml = ruamel.yaml.YAML()

# local imports
from silicate_classes import PlanetDef, MantleDef, VolcanoDef
from atmosphere_class import AtmosphereDef

def param_set(obj, params):
    """
    Sets the parameters of a class instance obj using params presented as a dict.

    Sets the attributes of a class according to parameters set in a yaml file.

    Args:
        obj (class): the class being set up
        params (dict): the parameters the class is being updated with.

    Returns:
        None
    """

    def Test_Params(par):
        """
        Checks the param is present in the class.
        
        Checks against properties of the class so the input file parameter has to match
        the class property name.
        
        Args:
            par (str): A parameter from the params dict.
        
        Returns:
            True if the parameter is present in the class.
        """
        
        for item in inspect.getmembers(obj):
            if item[0] == par:
                return True
        return 0
    
    # force mpfr if set
    for par, val in params.items():
        if Test_Params(par):
        #     if isinstance(getattr(obj, par), type(gp.mpfr(0))):
        #         val = str(val)
        #         setattr(obj, par, gp.mpfr(val))
        #     else:
        #         setattr(obj, par, val)
            setattr(obj, par, val)
        else:
            sys.exit("Warning: %s not a valid calculation parameter" % par)

def readin_planet_mantle(p, m):
    """
    Creates & sets up instances of the Planet & Mantle classes.

    Sets up instances of the two classes, checks the chemistry of the mantle against
    the list of allowed names (the ones we have data for).
    If a species is present in the mantle chemistry EVolve does not recognise, 
    the programme exits.

    Args:
        p (yaml): opened planet.yaml file
        m (yaml): opened mantle.yaml file

    Returns:
        planet (class): instance of the PlanetDef class
        mantle (class): instance of the MantleDef class

    """

    # creates a dictionary of the planet parameters from planet.yaml
    p_dat = yaml.full_load(p)

    # setup run definitions
    planet = PlanetDef()
    param_set(planet, p_dat)  # sets planetary parameters from the input file in the planet instance.
    
    # setup mantle definitions
    mantle = MantleDef()    
    m_dat = yaml.full_load(m)
    
    # Checks that only oxides we have data for are being listed as part of the melt.
    allowed_names = ['sio2', 'tio2', 'al2o3', 'feo', 'fe2o3', 'mno',\
                     'mgo', 'cao', 'na2o', 'k2o', 'p2o5', 'nio', 'cr2o3']
    
    for chem in m_dat['composition'].keys():
        if chem.lower() in allowed_names:
            pass
        else:
            exit(f"{chem} does not match an element in the allowed list (case insensitive):\n{allowed_names}")
    
    param_set(mantle, m_dat)
    mantle.init_store()

    return planet, mantle

#  -------------------------------------------------------------------------
def readin_volc(f, atm):
    """
    Creates & sets up an instance of the Volcano class.

    Args:
        f (dict): locations of the input files for EVo from input.yaml.
        atm (class): instance of the AtmosphereDef class

    Returns:
        volc (class): instance of the VolcanoDef class

    """
    
    volc = VolcanoDef(atm)
    param_set(volc, f)

    return volc

#  -------------------------------------------------------------------------
def readin_atm(f, planet, mantle, fchem):
    """
    Creates & sets up instance of the AtmosphereDef class.

    Sets up instance of theclass, updates with data from atmosphere.yaml then creates
    the atmosphere_out.csv file filled with the initial atmospheric conditions.
    Edits the config files for FastChem to match those in the inputs folder.

    Args:
        f (yaml): opened atm.yaml file

    Returns:
        atm (class): instance of the AtmosphereDef class

    """

    data = yaml.full_load(f)
    
    atm = AtmosphereDef()
    param_set(atm, data)
    fastchem = atm.initial_comp(planet, mantle, fchem)
    # edit config and parameter files for fastchem to those set in the inputs file of evolve
    atm.fc_setup()

    return atm, fastchem

#  -------------------------------------------------------------------------
def readin(input_data, fchem=False, nocrust=False):
    """
    Runs all the readin scripts.

    Takes the input data and seperates into relevent model sections.
    Creates instances of each class and updates them with info from
    these input files.

    Args:
        input_data (yaml): YAML file with the locations of all the input files listed.
        fchem (bool): If True, EVo output will be run through FastChem.

    Returns:
        Single instance each of PlanetDef, MantleDef, AtmosphereDef + VolcanoDef.
    """
    
    "Takes file locations from inputs.yaml and splits to allow setup from each file/location."
    
    with open(input_data, 'r') as input_files:
        data = yaml.full_load(input_files)

    p_dat= data['evolve']['planet']
    m_dat = data['evolve']['mantle']
    a_dat = data['evolve']['atmosphere']

    print("Reading in from:")
    print("Planet file:", p_dat)
    print("Mantle file:", m_dat)
    print("Atmosphere file:", a_dat,"\n")

    # Delete any existing output files
    out_fs = ['mantle_out.csv', 'volc_out.csv', 'fc_out.csv', 'fc_out_h2loss.csv', 'volc_mf_out.csv', 'atmosphere_out.csv', 'fc_shelloutput.txt', 'fc_input.csv']
    for f in out_fs:
        if os.path.isfile(f'outputs/{f}'):
            os.remove(f'outputs/{f}')
    
    # Save the input files to outputs/settings.txt
    writeout_settings(p_dat, m_dat, a_dat, fchem, nocrust)
    
    # planet & mantle
    with open(p_dat,"r") as p:
        with open(m_dat, "r") as m:
            planet, mantle = readin_planet_mantle(p, m)

    # atmosphere, if initial atmosphere is given
    if a_dat != 'None':
        with open(a_dat,"r") as a:
            atm, fchem = readin_atm(a, planet, mantle, fchem)
    else:
        atm = None
    
    if fchem == True:
        atm.fastchem == True

    # volcanic system
    evo_dat = data['evo']
    volc = readin_volc(evo_dat, atm)

    return planet, mantle, atm, volc, fchem

def writeout_settings(p, m, a, fchem, nocrust):
    """
    Saves all the input settings to a text file.

    The planet, mantle and atm.yaml files are saved to settings.txt so when
    outputs are saved the initial settings can be stored.
    """
    
    with open(p,"r") as planet:
        planet_file = ryaml.load(planet)
        
    with open(m,"r") as mantle:
        mantle_file = ryaml.load(mantle)
        
    with open(a,"r") as atm:
        atm_file = ryaml.load(atm)

    with open('inputs/env.yaml',"r") as env:
        env_file = ryaml.load(env)
    
    if not os.path.exists('outputs'):
        os.makedirs('outputs')
    
    with open('outputs/settings.txt', "w+") as settings_f:
        settings_f.write('OPTIONS:\n')
        settings_f.write(f'FastChem: {fchem}\n')
        settings_f.write(f'No crust: {nocrust}\n')
        settings_f.write(f'H2 escape: None available in this release')
        settings_f.write("\n")
        settings_f.write("\n")
        settings_f.write('# PLANET.YAML SETTINGS:')
        settings_f.write("\n")
        ryaml.dump(planet_file, settings_f)
        settings_f.write("# "+"-"*60)
        settings_f.write("\n")        
        settings_f.write("\n")
        settings_f.write('# MANTLE.YAML SETTINGS:')
        settings_f.write("\n")
        ryaml.dump(mantle_file, settings_f)
        settings_f.write("# "+"-"*60)
        settings_f.write("\n") 
        settings_f.write("\n")
        settings_f.write('# ATM.YAML SETTINGS:')
        settings_f.write("\n")
        ryaml.dump(atm_file, settings_f)
        settings_f.write("# "+"-"*60)
        settings_f.write("\n") 
        settings_f.write("\n")
        settings_f.write('# ENV.YAML SOLUBILITY MODELS:')
        settings_f.write("\n")
        settings_f.write(f"H2O_MODEL: {env_file['H2O_MODEL']}\n")
        settings_f.write(f"H2_MODEL: {env_file['H2_MODEL']}\n")
        settings_f.write(f"C_MODEL: {env_file['C_MODEL']}\n")
        settings_f.write(f"CO_MODEL: {env_file['CO_MODEL']}\n")
        settings_f.write(f"CH4_MODEL: {env_file['CH4_MODEL']}\n")
        settings_f.write(f"SULPHIDE_MODEL: {env_file['SULFIDE_CAPACITY']}\n")
        settings_f.write(f"SULPHATE_MODEL: {env_file['SULFATE_CAPACITY']}\n")
        settings_f.write(f"N_MODEL: {env_file['N_MODEL']}\n")
        settings_f.write("# "+"-"*60)
