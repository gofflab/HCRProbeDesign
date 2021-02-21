'''
@authors:
    Marshall J. Levesque
    Arjun Raj
    Daniel Wei
'''

import string
import math
import re
import array

def containsAny(astring, aset):
    # Check whether 'str' contains ANY of the chars in 'set'
    # http://code.activestate.com/recipes/65441-checking-whether-a-string-contains-a-set-of-chars/
    return 1 in [c in astring for c in aset]

def gibbs(dH,dS,temp=37):
    """ Calc Gibbs Free Energy in cal/mol from enthaply, entropy, and temperature

    Arguments:
    dH -- enthalpy in kcal/mol
    dS -- entropy in cal/(mol * Kelvin)
    temp -- temperature in celcius (default 37 degrees C)
    """
    return dH*1000 - (temp+273.15)*dS  # cal/mol


def init_rna_dna():
    """Return [enthalpy, entropy] list in kcal/mol and cal/(mol*Kelvin) for RNA/DNA
    duplex initiation. Values from Sugimoto et al 1995
    """
    initH = 1.9  # kcal/mol
    initS = -3.9 # cal/(mol * Kelvin)
    return [initH, initS]


def stacks_rna_dna(inseq):
    """Calculate RNA/DNA base stack thermodynamic values (Sugimoto et al 1995)

    Sugimoto 95 parameters for RNA/DNA Hybridization (Table 3)
    "Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid
    Duplexes" in Biochemistry 1995

    Input Arguments:
    inseq -- RNA sequence of the RNA/DNA hybrid ( 5'->3' uracil->thymidine)

    Return [enthalpy, entropy] list in kcal/mol and cal/(mol*Kelvin)
    """

    # uracil->thymidine
    delH = {'aa':-7.8, 'ac':-5.9, 'ag':-9.1, 'at':-8.3,
            'ca':-9.0, 'cc':-9.3, 'cg':-16.3,'ct':-7.0,
            'ga':-5.5, 'gc':-8.0, 'gg':-12.8,'gt':-7.8,
            'ta':-7.8, 'tc':-8.6, 'tg':-10.4,'tt':-11.5}  # kcal/mol

    delS = {'aa':-21.9, 'ac':-12.3, 'ag':-23.5, 'at':-23.9,
            'ca':-26.1, 'cc':-23.2, 'cg':-47.1, 'ct':-19.7,
            'ga':-13.5, 'gc':-17.1, 'gg':-31.9, 'gt':-21.6,
            'ta':-23.2, 'tc':-22.9, 'tg':-28.4, 'tt':-36.4}  # cal/(mol*Kelvin)

    # sum enthalpy and entropy of RNA-DNA base stacks
    dH = sum([delH[inseq[i:i+2]] for i in range(len(inseq)-1)]) # kcal/mol
    dS = sum([delS[inseq[i:i+2]] for i in range(len(inseq)-1)]) # cal/(mol*Kelvin)

    return [dH, dS]

def init_dna_dna(inseq):
    """Return [enthalpy, entropy] list with units kcal/mol and cal/(mol*Kelvin)
    for DNA/DNA duplex initiation for the input DNA sequence (actg 5'->3').
    Values from SantaLucia 1998. Argument is DNA
    """
    initH = 0  # kcal/mol
    initS = 0  # cal/(mol*Kelvin)

    if (inseq[0] == 'c') or (inseq[0] == 'g'):
        initH += 0.1
        initS += -2.8
    else:
        initH += 2.3
        initS += 4.1

    if (inseq[-1] == 'c') or (inseq[-1] == 'g'):
        initH += 0.1
        initS += -2.8
    else:
        initH += 2.3
        initS += 4.1

    return [initH, initS]


def stacks_dna_dna(inseq, temp=37):
    """Calculate thermodynamic values for DNA/DNA hybridization.

    Input Arguments:
    inseq -- the input DNA sequence of the DNA/DNA hybrid (5'->3')
    temp  -- in celcius for Gibbs free energy calc (default 37degC)
    salt  -- salt concentration in units of mol/L (default 0.33M)

    Return [enthalpy, entropy] list in kcal/mol and cal/(mol*Kelvin)
    """
    # SantaLucia 98 parameters for DNA Hybridization (Table 2)
    delH = {'aa':-7.9, 'ac':-8.4, 'ag':-7.8, 'at':-7.2,
            'ca':-8.5, 'cc':-8.0, 'cg':-10.6,'ct':-7.8,
            'ga':-8.2, 'gc':-9.8, 'gg':-8.0, 'gt':-8.4,
            'ta':-7.2, 'tc':-8.2, 'tg':-8.5, 'tt':-7.9}  # kcal/mol

    delS = {'aa':-22.2, 'ac':-22.4, 'ag':-21.0, 'at':-20.4,
            'ca':-22.7, 'cc':-19.9, 'cg':-27.2, 'ct':-21.0,
            'ga':-22.2, 'gc':-24.4, 'gg':-19.9, 'gt':-22.4,
            'ta':-21.3, 'tc':-22.2, 'tg':-22.7, 'tt':-22.2}  # cal/(mol*Kelvin)

    # sum enthalpy and entropy of DNA-DNA base stacks
    dH = sum([delH[inseq[i:i+2]] for i in range(len(inseq)-1)]) # kcal/mol
    dS = sum([delS[inseq[i:i+2]] for i in range(len(inseq)-1)]) # cal/(mol*Kelvin)

    return [dH, dS]


def salt_adjust(delG,nbases,saltconc):
    """Adjust Gibbs Free Energy from 1M Na+ for another concentration

    Arguments:
    delG -- Gibbs free energy in kcal/mol
    nbases --  number of bases in the sequence
    saltconc -- desired Na+ concentration for new Gibbs free energy calculation

    Equation 7 SantaLucia 1998
    """
    return delG - 0.114*nbases*math.log(saltconc)


def overhang_rna(inseq,end):
    """Return Gibbs free energy at 37degC (in kcal/mol) contribution from single
    base overhang in RNA/RNA duplex.

    Arguments:
    inseq - 2bp RNA sequence (5' -> 3') uracil->thymidine
    end - specifies which end the over hang is on (valid values: 3 or 5)

    Table 3 in Freier et al, Biochemistry, 1986
    """
    # Free energy in kcal/mol for RNA/RNA 1M NaCl, 37 degrees celcius
    if (end == 5):
        dGoh = {'aa':-0.3, 'ac':-0.5, 'ag':-0.2, 'at':-0.3,
                'ca':-0.3, 'cc':-0.2, 'cg':-0.3, 'ct':-0.2,
                'ga':-0.4, 'gc':-0.2, 'gg':-0.0, 'gt':-0.2,
                'ta':-0.2, 'tc':-0.1, 'tg':-0.0, 'tt':-0.2}
    elif (end == 3):
        dGoh = {'aa':-0.8, 'ac':-0.5, 'ag':-0.8, 'at':-0.6,
                'ca':-1.7, 'cc':-0.8, 'cg':-1.7, 'ct':-1.2,
                'ga':-1.1, 'gc':-0.4, 'gg':-1.3, 'gt':-0.6,
                'ta':-0.7, 'tc':-0.1, 'tg':-0.7, 'tt':-0.1}

    return dGoh[inseq]


def overhang_dna(inseq,end):
    """Return Gibbs free energy at 37degC (in kcal/mol) contribution from single
    base overhang in DNA/DNA duplex.

    Arguments:
    inseq - 2bp DNA sequence (5' -> 3')
    end - specifies which end the over hang is on (valid values: 3 or 5)

    Table 2 in Bommarito, S. (2000). Nucleic Acids Research
    """

    # Free energy in kcal/mol for DNA/DNA 1M NaCl, 37 degrees celcius
    if (end == 5):
        dGoh = {'aa':-0.51, 'ac':-0.96, 'ag':-0.58, 'at':-0.50,
                'ca':-0.42, 'cc':-0.52, 'cg':-0.34, 'ct':-0.02,
                'ga':-0.62, 'gc':-0.72, 'gg':-0.56, 'gt':-0.48,
                'ta':-0.71, 'tc':-0.58, 'tg':-0.61, 'tt':-0.10}
    elif (end == 3):
        dGoh = {'aa':-0.12, 'ac':+0.28, 'ag':-0.01, 'at':+0.13,
                'ca':-0.82, 'cc':-0.31, 'cg':-0.01, 'ct':-0.52,
                'ga':-0.92, 'gc':-0.23, 'gg':-0.44, 'gt':-0.35,
                'ta':-0.48, 'tc':-0.19, 'tg':-0.50, 'tt':-0.29}

    return dGoh[inseq]


def melting_temp(dH,dS,ca,cb,salt):
    # calculation for total concentration of nucleic acid for non-self-complementary pairs
    ct = ca - cb/2 # SantaLucia 98

    # calculation for melting temperature - SantaLucia 98 Equation 3
    tm = dH*1000/(dS + (1.9872 * math.log(ct))) + (16.6 * math.log10(salt)) - 273.15
    return tm
