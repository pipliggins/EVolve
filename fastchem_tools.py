import numpy as np
import os, csv, ctypes, io, tempfile, sys
from contextlib import contextmanager

from fastchem.python import pyfastchem
import conversions as cnvs

def saveChemistryOutput(filename, temperature, pressure, fastchem_output, fastchem, tstep, tracked_species=None):
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
  
    #total gas particle number density from the ideal gas law 
    #used to convert the number densities to mixing ratios
    const_k = 1.3806504e-16     # boltzmann's constant in erg K-1, to match exactly the C++ version. Internet & astropy say 1.380649e-16
    
    gas_number_density = pressure*1e6 / (const_k * temperature)

    if tracked_species != None:
        tracked_sp = {cnvs.evo2fc_names(sp):j for sp, j in zip(tracked_species, range(len(tracked_species)))}
    else:
        tracked_sp = None

    #calculate the mixing ratios from the number densities
    mixing_ratios = np.array(fastchem_output.number_densities) / gas_number_density

    fc_start = [tstep, pressure, temperature, fastchem_output.total_element_density[0], gas_number_density, fastchem_output.mean_molecular_weight[0]]
    
    file_exists = os.path.isfile('outputs/' + filename)

    if file_exists:
        with open('outputs/'+filename, 'a') as f_out:
            write_fout = csv.writer(f_out, delimiter='\t')
            write_fout.writerow(fc_start + list(mixing_ratios[0]))
    else:
        titles = ['t_step', 'P (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (g/mol)'] + [fastchem.getSpeciesSymbol(j) for j in range(fastchem.getSpeciesNumber())]
        
        with open('outputs/'+filename, 'w') as f_out:
            write_fout = csv.writer(f_out, delimiter='\t')
            write_fout.writerow(titles)
            write_fout.writerow(fc_start + list(mixing_ratios[0]))

    if tracked_species != None:
        for j in range(fastchem.getSpeciesNumber()):
            if fastchem.getSpeciesSymbol(j) in tracked_sp.keys():
                tracked_sp[fastchem.getSpeciesSymbol(j)] = mixing_ratios[0, j]

    for evo_n in tracked_species:
        tracked_sp[evo_n] = tracked_sp.pop(cnvs.evo2fc_names(evo_n))
    
    return fastchem_output.mean_molecular_weight[0], tracked_sp

libc = ctypes.CDLL(None)
c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')

@contextmanager
def stdout_redirector(stream):
    """
    NO LONGER USED: Really messed around with printing in the rest of the programme, so edited the base C instead.
    
    Redirects the terminal output to the variable 'stream'.

    Used to stop the 'neglected' terms being printed by fastchem when
    using the python module. Taken from
    https://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/
    
    Usage:
        f = io.BytesIO()
        with stdout_redirector(f):
            # run commands without any terminal output printing.
    """
    
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stdout_fd)
