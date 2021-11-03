"""
Stores universally useful functions such as unit conversions.

Contains functions to:
    Normalise a composition
    Convert a composition from weight fractions to mole fractions
    Convert standard species names to those used in FastChem outputs
    Perform standard mathematical rounding.
    
"""

import math
import constants as cnst
import constants as cnst

def norm(c, f=1., cons=None):
    """
    Normalises the input composition to value defined by f (or 1)

    Sum to f. If constants are present (cons), then these are held to a constant
    value and the rest are normalised.

    Args:
        c (dict): Dictionary of chemical species and their abundances
        f (float): Number to normalise to.
        cons (lst of sts): List of species to hold constant.

    Returns:
        tmp (dict): Normalised version of c

    """
    
    sm = 0
    tmp = {}
    constant = 0

    if cons:
        for x in c:
            if x in cons:
                constant += c[x]
            else:
                sm += c[x]
    else:
        for x in c:
            sm += c[x]

    if cons:
        for x in c:
            if x in cons:
                tmp[x] = c[x]
            else:
                tmp[x] = (f-constant) * c[x] / sm
    else:
        for x in c:
            tmp[x] = f * c[x] / sm

    return tmp

def wt2mol(c, *args):
    """
    Convert weight fractions to mol fractions.

    Returns a composition dict converted from weight to mole fractions.
    If args, only the species listed in args are returned.

    Args:
        c (dict): Composition as a weight fraction
        args (str): Specific species wanted to return.

    Returns:
        mol (dict): Composition as mol fractions.

    """
    
    mol = {}
    sm = 0.0

    for x in c:
        if x == 'o(fe)':
            mol[x] = c[x] / cnst.m['o']
            sm += c[x] / cnst.m['o']
        else:
            mol[x] = c[x] / cnst.m[x.lower()]
            sm += c[x] / cnst.m[x.lower()]

    for x in c:
        mol[x] = mol[x] / sm

    if args:  # should be able to specify one element from the list to return if necessary.
        for arg in args:
            return mol[arg]
    else:
        return mol

def mol2wt(c, i=-1):
    """
    Convert a composition from mol fractions to weight fractions.

    Args:
        c (dict): Composition as a mole fractions
        i (int): Index specifier if c is a dict made up of lists

    Returns:
        cw (dict): re-normalised composition as weight fractions.
    """

    cw = {}
    sm = 0

    for x in c:
        try:
            cw[x] = c[x][i] * cnst.m[x.lower()]
            sm += c[x][i] * cnst.m[x.lower()]
        except:
            cw[x] = c[x] * cnst.m[x.lower()]
            sm += c[x] * cnst.m[x.lower()]

    for x in cw:
        cw[x] = cw[x] / sm

    return cw

def mass2wt_frac(c):
    """
    Convert a composition from absolute masses to weight fractions.

    Args:
        c (dict): composition of a phase as a set of absolute masses.

    Returns:
        cw (dict): Composition as weight fractions.
    """

    cw = {}
    sm = 0

    for s in c:
        sm += c[s]

    for s in c:
        cw[s] = c[s]/sm

    return cw

def wt_frac2mass(c, mass):
    """
    Convert a composition to actual masses.

    Args:
        c (dict): composition of a phase as weight fractions.
        mass (float): Mass of the phase who's composition is in c.

    Returns:
        cmass (dict): Composition as absolute masses.
    """

    cmass = {}

    for s in c:
        cmass[s] = c[s]*mass
    
    return cmass

def evo2fc_names(name):
    """
    Converts the chemical names as used in EVo and EVolve to those in FastChem

    e.g: Converts HCN -> C1H1N1

    Args:
        name (str): Chemical species name in standard notation.

    Returns:
        with_nums (str): Converted name as a string.

    """

    letters = list(name)  # HCN -> H, C, N

    length = len(letters)

    i = 0
    while i <= length:
        if letters[i] == letters[-1]:   # if the last letter
            if not letters[-1].isnumeric() and not (letters[-1].islower() and len(letters)==2):
                letters.append('1')
                break
            else:
                break
        elif letters[i].isnumeric() or letters[i+1].isnumeric():
            i += 1
            continue
        elif letters[i+1].islower():
            i += 1
            continue
        else:
            letters.insert(i+1, '1')
            length += 1  # length has now increased
            i += 1
    
    with_nums = ''.join(letters)    # H1C1N1

    if len(with_nums) > 2:
        n=2  # size of each string, so each sring is a species-number pair
        split_2 = [with_nums[i:i+n] for i in range(0, len(with_nums), n)]   # splits into 'H1', 'C1', 'N1'

        alphab = sorted(split_2)

        # return ''.join(alphab)    # C1H1N1
        final =  ''.join(alphab)    # C1H1N1

        if final != 'C1H1N1':
            return final
        else:
            return final+'_1'
    
    else:
        return with_nums
    
def round_proper(n, decimals=0):
    """
    Perfoms 'traditional' rounding to the specified number of d.p.

    Args:
        n (float): The number to be rounded.
        decimals (int): The number of decimal places to round to.

    Returns:
        n rounded to the correct number of decimal places.
    """
    
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier

def mf2fi(mf):
    """
    Convert mole fraction to f_i, where f_i = ni/na, and na=(nt-ni)
    
    Args:
        mf (float): mole fraction of gas in the atmosphere
    Returns:
        fi (float): the mixing ratio
    """
    
    return mf/(1-mf)    

def fmq2fo2(dfmq,t,p):
    """
    Convert an fo2 value relative to the FMQ buffer to log10(fo2).

    Args:
        dfmq (float): fo2 relative to fmq
        t (float): temperature (K)
        p (float): pressure (bar)
    
    Returns:
        fo2 (float): log10(fo2)
    """
    
    return(dfmq + (-25096.3/t + 8.735 + 0.11*(p-1)/t))

def iw2fo2(FO2, t, p):
    """
    Convert an fo2 value relative to the IW buffer to log10(fo2).

    Args:
        dfmq (float): fo2 relative to iw
        t (float): temperature (K)
        p (float): pressure (bar)
    
    Returns:
        fo2 (float): log10(fo2)
    """
    return(FO2 + (-27489/t + 6.702 + 0.055*(p-1)/t))

def nno2fo2(FO2, t, p):
    """
    Convert an fo2 value relative to the NiNiO buffer to log10(fo2).

    Args:
        dfmq (float): fo2 relative to nno
        t (float): temperature (K)
        p (float): pressure (bar)
    
    Returns:
        fo2 (float): log10(fo2)
    """
    return(FO2 + (-24930/t + 9.36 + 0.046*(p-1)/t))


def generate_fo2(dfo2, buffer, P, T):
    """
    Takes an fo2 value relative to any buffer, and returns oxygen fugacity.
    
    Args:
        dfo2 (float): fo2 relative to a rock buffer
        buffer (string): Rock buffer abbreviation (i.e. FMQ, NNO, IW)
        P (float): pressure (bar)
        T (float): temperature (K)

    Returns:
        fo2 (float): standard oxygen fugacity
    """

    if buffer == 'FMQ':
        fo2 = 10 ** fmq2fo2(dfo2, T, P)

    elif buffer == 'IW':
        fo2 = 10 ** iw2fo2(dfo2, T, P,)
 
    elif buffer == 'NNO':
        fo2 = 10 ** nno2fo2(dfo2, T, P)

    return fo2

def atomic_mass2mantle_volatiles(atomic_fracs):
        """
        Takes melt volatile contents as ppm H2O, CO2 etc and calculates the atomic
        weight fractions of H, S, C and N.

        Args:
            atomic_fracs (list): Containing the atomic weight fractions (H, C, S, N) in the melt as ppm

        Returns:
            tuple: (H2O, CO2, S, N) as mass fractions (ppm).
        """

        h, c, s, n = atomic_fracs

        h2o = (h/cnst.m['h'])*cnst.m['h2o']*0.5     # 2 moles of h for every one of H2O
        co2 = (c/cnst.m['c'])*cnst.m['co2']

        return h2o, co2, s, n

def calc_mmw(comp, i=-1):
    """
    Calculates the mean moleuclar weight of the gas phase

    Args:
        comp (dict): The composition pf the atmosphere as mol fractions.
        i (int): the index to use if comp is a dict of lists.

    Returns:
        sm (float): the mean molecular weight of the gas.
    """

    sm = 0

    for sp in comp.keys():
        try:
            sm += comp[sp]*cnst.m[sp.lower()]
        except:
            sm += comp[sp][i]*cnst.m[sp.lower()]

    return sm