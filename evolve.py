"""
The main module for EVolve.

A model of atmospheric buildup from a volcanic source, with inbuilt magma-gas chemistry from EVo, atmospheric chemistry from FASTCHEM and atmospheric losses.
EVolve takes a planet with an initial mantle chemistry, volatile content, and atmosphere and models the evolution of mantle and atmosphere as volcanic degassing moves volatiles
from the interior to the atmosphere.


Planet-wide settings such as mass, radius and stellar environment are set in planet.yaml
Mantle-specific settings such as the major oxide geochemistry, total volatile content and melt partitioning are set in mantle.yaml
Atmospheric constraints such as the initial size/thickness, mean molecular mass and gas chemistry are stored in atm.yaml

Up to 4 results files will be generated, stored in the 'outputs' file:
    mantle_out.csv - storing the fo2 and volatile content of the mantle for each timestep
    volc_out.csv - storing the final pressure step from EVo for each timestep, in EVo's units. 
    fc_out.csv - present if the fastchem option is selected, will store the volcanic input to the atmosphere
                         at each timestep after re-equilibration to atmospheric temperature and processing by fastchem, atmospheric constituents present as volume mixing ratios (actually mole fractions)
    volc_mf_out.csv - present if the fastchem option is NOT selected. Stores the volcanic input to the atmosphere from EVo, as mol fractions
    atmosphere_out.csv - stores the atmospheric pressure and composition at the surface at the end of each timestep, once the volcanic input has been mixed in. Composition in VMR.
                            Only species listed as being tracked will be stored here.

"""

#------------------------------------------------------------------------
# IMPORTS
#------------------------------------------------------------------------

# python main [ ensure these are on your system ]
import argparse

# bundled scripts
import readin
import messages as msgs

#------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------

def main(input_f, fastchem, nocrust):
    """
    The main function which runs EVolve.

    Args:
        input_f (string): the name of the inputs file which contains the locations of all the input files.
        fastchem (bool): If true, FastChem equilibrium chemistry is run.
        nocrust (bool): If true, all non-erupted volatiles and all melt is returned to the mantle 
                        immediately, and no crust is formed.
    """ 
    
    # Reads in settings from files
    planet, mantle, atm, volc, fastchem = readin.readin(input_f, fchem=fastchem, nocrust=nocrust)
    planet.setup_conservation(mantle, atm)

    t_step = 1

    while t_step <= planet.max_timestep:
    
        # generate melt
        volc.melt_volatiles, volc.melt_mass = mantle.generate_melt(planet, atm, t_step, nocrust=nocrust)

        # edit evo input files according to melt volatile load, fO2 etc
        mantle.set_evo_env(volc, atm, set_atomic_weights=planet.set_melt_atomic_weights)

        # Run evo with edited env input file    
        volc.run_evo(t_step)
        mantle.update_crust(planet, volc, atm, volc.melt_fracs, volc.WgT, volc.melt_mass, volc.melt_volatiles, t_step, nocrust=nocrust)
        atm.update_store(volc, volc.melt_mass)

        new_psurf, p_volc = atm.update_psurf(planet, volc)

        if fastchem == True:   
            # run the mixed atmosphere (volcanic emissions + pre-existing atmosphere) through fastchem.
            
            atm.p_surf = new_psurf
            atm.fc_cheminput = atm.unitconvert_atomicmass2fc()
            
            atm.run_fastchem(planet, t_step)
            atm.fc_saveatmosphere(planet, t_step)
        
        else:
            # add the evo result directly to the existing composition @ magmatic temperatures
            atm.atmos_input = volc.mol_fracs       
            new_comp, atm.p_surf = atm.mix_atmosphere(new_psurf, p_volc)
            
            atm.saveatmosphere(planet, new_comp, t_step)
        
        planet.check_conservation(mantle, atm, nocrust=nocrust)
        
        # Check for atmospheric chem convergence, or if one loop has been specified
        if planet.snapshot == True:
            msgs.exiting('Single run completed.')
        
        elif atm.check_evolution() == False:
            with open('outputs/settings.txt', "a") as settings_f:
                settings_f.write("\n")
                settings_f.write(f"INFO: Atmosphere appears to have reached steady state. Continuing...")
                settings_f.write("\n")
            
            t_step += 1
        
        else:
            t_step += 1
    
    msgs.exiting('Max number of iterations reached without atmosphere chemistry reaching steady state. EVolve has finished.')


if __name__ == "__main__":
       
    # Create the parser
    my_parser = argparse.ArgumentParser(prog='evolve', description='Run EVolve: build a volcanic atmosphere')

    # Add the arguments
    my_parser.add_argument('input_files',
                        metavar='inputs.yaml',
                        help='lists locations of each input file for evolve, evo and fastchem')
    
    my_parser.add_argument('--fastchem',
                        action='store_true',
                        help='Include fastchem equilibrium chemistry in the atmosphere')
    
    my_parser.add_argument('--nocrust',
                        action='store_true',
                        help='Return all unerupted volatiles and all melt volume to the mantle, no crust reservoir is formed.')

    # Parse in files
    args = my_parser.parse_args()

    input_f = args.input_files

    if args.fastchem == True:
        fastchem = True
    else:
        fastchem = False

    if args.nocrust == True:
        nocrust = True
    else:
        nocrust = False
    
    main(input_f, fastchem, nocrust)
