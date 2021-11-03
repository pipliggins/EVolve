import yaml
import ruamel.yaml
from itertools import product
from shutil import copytree
import argparse
from datetime import date

from evolve import main
ryaml = ruamel.yaml.YAML()

def amend_files(run, p, m, a):
    """
    Edits the EVolve input files.

    Args:
        run (dict): the current permutation which files should be updated to reflect
        p (dict): all the parameters to iterate through in the planet.yaml file
        m (dict): all the parameters to iterate through in the mantle.yaml file
        a (dict): all the parameters to iterate through in the atmosphere.yaml file
    
    Returns:
        None
    """

    for k, v in run.items():
        if k in p.keys():
            with open('inputs/planet.yaml', 'r') as f:
                p_file = ryaml.load(f)
                p_file[k] = v
                
            with open('inputs/planet.yaml', "w") as fw:
                ryaml.dump(p_file, fw)
            
            continue
        
        elif k in m.keys() or ('mantle_volatiles' in m.keys() and k in m['mantle_volatiles'].keys()) or ('composition' in m.keys() and k in m['composition'].keys()):
            if 'mantle_volatiles' in m.keys() and k in m['mantle_volatiles'].keys():
                with open('inputs/mantle.yaml', 'r') as f:
                    m_file = ryaml.load(f)
                    m_file['mantle_volatiles'][k] = v

            elif 'composition' in m.keys() and k in m['composition'].keys():
                with open('inputs/mantle.yaml', 'r') as f:
                    m_file = ryaml.load(f)
                    m_file['composition'][k] = v
            
            else:
                with open('inputs/mantle.yaml', 'r') as f:
                    m_file = ryaml.load(f)
                    m_file[k] = v
            
            with open('inputs/mantle.yaml', "w") as fw:
                ryaml.dump(m_file, fw)
            
            continue
        
        else:
            if 'composition' in a.keys() and k in a['composition'].keys():
                with open('inputs/atm.yaml', 'r') as f:
                    a_file = ryaml.load(f)
                    a_file['composition'][k] = v

            elif 'tracked_species' in a.keys() and k in a['tracked_species'].keys():
                with open('inputs/atm.yaml', 'r') as f:
                    a_file = ryaml.load(f)
                    a_file['tracked_species'][k] = v
            
            else:
                with open('inputs/atm.yaml', 'r') as f:
                    a_file = ryaml.load(f)
                    a_file[k] = v

            with open('inputs/atm.yaml', "w") as fw:
                ryaml.dump(a_file, fw)
                

def multirun(run_name, file, fastchem, nocrust):
    """
    The main function when parameter sweeping using EVolve.

    Finds all possible permutations of the parameters to be varied,
    specified in multirun.yaml, then runs each permutation by
    updating the input files then running evolve.

    All results from the parameter sweep are stored in a folder within 'results'.

    Args:
        run_name (str): The folder name for this parameter sweep, set on the commandline
        file (str): The yaml filename containing the parameters to sweep over
        fastchem (bool): If True, fastchem equilibrium chemistry is used in evolve.

    Returns:
        None
    """
    
    with open(file, 'r') as f:
        multirun_doc = yaml.full_load(f)

        planet = multirun_doc['planet']
        mantle = multirun_doc['mantle']
        atm = multirun_doc['atmosphere']

    # then loop, editing the input files and running evolve.

    merged_dicts = {}
    for group in [planet, mantle, atm]:
        if group != None:
            new_group = {}
            for k, v in group.items():
                if isinstance(v, dict): # any nested dictionares are included in the dict for calculating permutations
                    for a, b in v.items():
                        new_group[a] = b  
                else:
                    new_group[k] = v                   
            merged_dicts.update(new_group)

    keys, values = zip(*merged_dicts.items())
    permutations_dicts = [dict(zip(keys, v)) for v in product(*values)]

    for run in permutations_dicts:
                
        amend_files(run, planet, mantle, atm)

        vals = list(run.values())
        iteration_name = ''
        for i in vals:
           iteration_name += str(i)+'_' 

        # Runs EVolve
        try:
            main('inputs.yaml', fastchem, nocrust)
            
            # Copies the output folder into a separate folder within 'results'
            copytree('outputs', f'results/{run_name}/{iteration_name}')

        except:
            copytree('outputs', f'results/{run_name}/{iteration_name}')
            continue        


if __name__ == "__main__":
       
    # Create the parser
    my_parser = argparse.ArgumentParser(prog='evolve_multirun', description='Run a set of EVolve runs')

    # Add the arguments
    my_parser.add_argument('run_name',
                        metavar=f'multirun_{date.today().strftime("%d-%m-%y")}',
                        help='name of the run set to label the results file')
    
    my_parser.add_argument('multirun_file',
                        metavar='multirun.yaml',
                        help='contains the properties to iterate over')
    
    my_parser.add_argument('--fastchem',
                        action='store_true',
                        help='Include fastchem equilibrium chemistry in the atmosphere')

    my_parser.add_argument('--nocrust',
                        action='store_true',
                        help='Return all unerupted volatiles and all melt volume to the mantle, no crust reservoir is formed.')

    # Parse in files
    args = my_parser.parse_args()

    folder_name = args.run_name
    multirun_f = args.multirun_file

    if args.fastchem == True:
        fastchem = True
    else:
        fastchem = False

    if args.nocrust == True:
        nocrust = True
    else:
        nocrust = False
   
    multirun(folder_name, multirun_f, fastchem, nocrust)