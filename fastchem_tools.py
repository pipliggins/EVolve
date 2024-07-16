import numpy as np
import os, csv, ctypes, io, tempfile, sys
from contextlib import contextmanager

import conversions as cnvs


def saveChemistryOutput(
    filename,
    temperature,
    pressure,
    fastchem_output,
    fastchem,
    tstep,
    tracked_species=None,
):
    """
    Mimics the file save process from FastChem's C++ library.

    Edited to return the mean molecular weight, and the mole fractions of species
    being tracked in evolve.

    Args:
        filename (str): file name/path being written out to.
        temperature (float): atmospheric temperature
        pressure (float): surface pressure
        fastchem_output (class): fastchem output data from python run
        fastchem (class): fastchem python object
        tstep (float): the timestep in evolve
        tracked_species (lst): the species being tracked in evolve.

    Returns:
        fastchem_output.mean_molecular_weight (float): mean moleuclar weight
        tracked_sp (dict): tracked species with associated mole fractions.
    """

    # total gas particle number density from the ideal gas law
    # used to convert the number densities to mixing ratios
    const_k = 1.3806504e-16  # boltzmann's constant in erg K-1, to match exactly the C++ version. Internet & astropy say 1.380649e-16

    gas_number_density = pressure * 1e6 / (const_k * temperature)

    if tracked_species != None:
        tracked_sp = {
            cnvs.evo2fc_names(sp): j
            for sp, j in zip(tracked_species, range(len(tracked_species)))
        }
    else:
        tracked_sp = None

    # calculate the mixing ratios from the number densities
    mixing_ratios = np.array(fastchem_output.number_densities) / gas_number_density

    fc_start = [
        tstep,
        pressure,
        temperature,
        fastchem_output.total_element_density[0],
        gas_number_density,
        fastchem_output.mean_molecular_weight[0],
    ]

    file_exists = os.path.isfile("outputs/" + filename)

    if file_exists:
        with open("outputs/" + filename, "a") as f_out:
            write_fout = csv.writer(f_out, delimiter="\t")
            write_fout.writerow(fc_start + list(mixing_ratios[0]))
    else:
        titles = [
            "t_step",
            "P (bar)",
            "T (K)",
            "n_<tot> (cm-3)",
            "n_g (cm-3)",
            "m (g/mol)",
        ] + [
            fastchem.getGasSpeciesSymbol(j)
            for j in range(fastchem.getGasSpeciesNumber())
        ]

        with open("outputs/" + filename, "w") as f_out:
            write_fout = csv.writer(f_out, delimiter="\t")
            write_fout.writerow(titles)
            write_fout.writerow(fc_start + list(mixing_ratios[0]))

    if tracked_species != None:
        for j in range(fastchem.getGasSpeciesNumber()):
            if fastchem.getGasSpeciesSymbol(j) in tracked_sp.keys():
                tracked_sp[fastchem.getGasSpeciesSymbol(j)] = mixing_ratios[0, j]

    for evo_n in tracked_species:
        tracked_sp[evo_n] = tracked_sp.pop(cnvs.evo2fc_names(evo_n))

    return fastchem_output.mean_molecular_weight[0], tracked_sp
