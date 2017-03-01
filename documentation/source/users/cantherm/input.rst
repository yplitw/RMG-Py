*************************************************************************************
Creating Input Files for Thermodynamics and High-Pressure Limit Kinetics Computations
*************************************************************************************

Syntax
======

The format of CanTherm input files is based on Python syntax. In fact, CanTherm
input files are valid Python source code, and this is used to facilitate 
reading of the file. 

Each section is made up of one or more function calls, where parameters are 
specified as text strings, numbers, or objects. Text strings must be wrapped in
either single or double quotes.

The following is a list of all the functions of a CanTherm input file:

=========================== =========================================================
Function                    Description
=========================== =========================================================
``modelChemistry``          Level of theory from quantum chemical calculations
``frequencyScaleFactor``    A factor by which to scale all frequencies
``useHinderedRotors``       ``True`` if hindered rotors are used, ``False`` if not
``useBondCorrections``      ``True`` if bond corrections are used, ``False`` if not
``species``                 Contains parameters for non-transition states
``transitionState``         Contains parameters for transition state(s)
``reaction``                Required for performing kinetic computations
``statmech``                Loads statistical mechanics parameters
``thermo``                  Performs a thermodynamics computation
``kinetics``                Performs a high-pressure limit kinetic computation
=========================== =========================================================

********* missing: network, pressureDependence ************

Note that ``frequencyScaleFactor`` can defined glabally or per ``species()`` or
``transitionState()`` function.


Model Chemistry
===============

The first item in the input file should be a ``modelChemistry()`` function,
which accepts a string describing the model chemistry.

CanTherm uses this information to adjust the computed energies to the usual gas-phase reference
states by applying atom, bond and spin-orbit coupling energy corrections. This is particularly
important for ``thermo()`` calculations (see below). Note that the user must specify under the
``species()`` function the type and number of atoms and bonds for CanTherm to apply these corrections.
The example below specifies CBS-QB3 as the model chemistry::

    modelChemistry("CBS-QB3")

Currently, atomization energy corrections (AEC), bond corrections (BC), and spin orbit correction (SOC) are available for the following model chemistries:

================================================ ===== ==== ====
Model Chemistry                                  AEC   BC   SOC     
================================================ ===== ==== ====
``'CBS-QB3'``                                     v    v    v
``'G3'``                                          v         v
``'M08SO/MG3S*'``                                 v         v
``'Klip_1'``                                      v         v
``'Klip_2'`` *uses QCI(tz,qz) values*             v         v
``'Klip_3'`` *uses QCI(dz,qz) values*             v         v
``'Klip_2_cc'`` *uses CCSD(T)(tz,qz) values*      v         v
``'CCSD-F12/cc-pVDZ-F12'``                        v         v
``'CCSD(T)-F12/cc-pVDZ-F12_H-TZ'``                v         v
``'CCSD(T)-F12/cc-pVDZ-F12_H-QZ'``                v         v
``'CCSD(T)-F12/cc-pVnZ-F12'``, *n = D,T,Q*        v    v    v
``'CCSD(T)-F12/cc-pVDZ-F12_noscale'``             v         v
``'CCSD(T)-F12/cc-pCVnZ-F12'``, *n = D,T,Q*       v         v
``'CCSD(T)-F12/aug-cc-pVnZ-F12'``, *n = D,T,Q*    v         v
``'B-CCSD(T)-F12/cc-pVnZ-F12'``, *n = D,T,Q*      v         v
``'B-CCSD(T)-F12/cc-pCVnZ-F12'``, *n = D,T,Q*     v         v
``'B-CCSD(T)-F12/aug-cc-pVnZ-F12'``, *n = D,T,Q*  v         v
``'DFT_G03_b3lyp'``                               v    v    v
``'DFT_ks_b3lyp'``                                          v
``'DFT_uks_b3lyp'``                                         v
``'G03_PBEPBE_6-311++g_d_p'``                     v         v
``'MP2_rmp2_pVnZ'``, *n = D,T,Q*                  v         v
``'FCI/cc-pVnZ'``, *n = D,T,Q*                    v         v
``'BMK/cbsb7'``                                   v    v    v
``'BMK/6-311G(2d,d,p)'``                          v    v    v
``'B3LYP/6-311+G(3df,2p)'``                            v
================================================ ===== ==== ====

Notes: In ``'M08SO/MG3S*'`` the grid size used in the [QChem] electronic structure calculation utilizes 75 radial points and 434 angular points. ``'DFT_G03_b3lyp'`` is a B3LYP calculation with a moderately large basis set.

Species
=======

Each species of interest must be specified using a ``species()`` function,
which accepts the following parameters:

====================== =========================================================
Parameter              Description
====================== =========================================================
``label``              A unique string label used as an identifier
``geometry``           The path to the quantum chemistry output file containing the optimized geometry
``frequencies``        The path to the quantum chemistry output file containing the computed frequencies
``externalSymmetry``   The external symmetry number for rotation
``linear``             ``True`` if the molecule is linear, ``False`` if not
``rotors``             A list of :class:`HinderedRotor()` objects describing the hindered rotors (optional)
``atoms``              Type and number of atoms in the species
``bonds``              Type and number of bonds in the species
====================== =========================================================

These parameters can be manually entered in the input file if known from the literture (e.g.,
from `CCCBDB <http://cccbdb.nist.gov/>`_), or automatically parsed by pointing to the output
file of a quantum chemistry calculation. In case of parsing, the externalSymmetry (see `thermo at CCCBDB <http://cccbdb.nist.gov/thermo.asp>`_ must still be manually entered.

Allowed atom symbols for the ``atoms`` parameter are 
``'C'``, ``'N'``, ``'O'``, ``'S'``, ``'P'``, and ``'H'``. For example, for formaldehyde we would write::

    atoms = {'C': 1, 'O': 1, H': 2}

Allowed bond types for the ``bonds`` parameter are, e.g., ``'C-H'``, ``'C-C'``, ``'C=C'``, ``'N-O'``, ``'C=S'``, ``'O=O'``, ``'C#N'``...

(The element order is according to the order specified under "allowed atom symbols" above, i.e., ``'C=N'`` is correct while ``'N=C'`` is incorrect. Use ``-``/``=``/``#`` to denote a single/double/triple bond, respectively). For example, for formaldehyde we would write::

    bonds = {'C=O': 1, 'C-H': 2}


Each :class:`HinderedRotor()` object requires the following parameters:

====================== =========================================================
Parameter              Description
====================== =========================================================
``scanLog``            The path to the Gaussian/Qchem log file containing the scan
``pivots``             The indices of the atoms in the hindered rotor torsional bond
``top``                The indices of all atoms on one side of the torsional bond (including the pivot atom)
``symmetry``           The symmetry number for the torsional rotation
====================== =========================================================

The following is an example of a typical species item, based on ethane::

    species(
        label = 'ethane',
        geometry = 'ethane_cbs.log',
        frequencies = 'ethane_cbs.log',
        extSymmetry = 2,
        linear = False,
        rotors = [
            HinderedRotor(scanLog='ethane_scan_1.log', pivots=[0,4], top=[0,1,2,3], symmetry=3),
        ]
        atoms = {'C': 2, 'H': 6},
        bonds = {'C-C': 1, 'C-H': 6},
    )

Note that the atoms identified within the rotor section should correspond to the indicated geometry. 

Transition State
================

Transition state(s) are only required when performimg kinetics computations.
Each transition state of interest must be specified using a ``transitionState()``
function, which is analogous to the ``species()`` function described above.

The following is an example of a typical transition state item::

    transitionState(
        label = 'TS1', 
        geometry = 'H+C2H4.log', 
        frequencies = 'H+C2H4.log', 
        extSymmetry = 2,
        linear = False, 
        rotors = [],
        atoms = {'C': 2, 'H': 5},
        bonds = {'C-C': 1, 'C-H': 5},
    )

Reaction
========

This is only required if you wish to perform a kinetics computation.
Each reaction of interest must be specified using a ``reaction()`` function,
which accepts the following parameters: 

====================== =========================================================
Parameter              Description
====================== =========================================================
``label``              A unique string label used as an identifier
``reactants``          A list of strings indicating the labels of the reactant species
``products``           A list of strings indicating the labels of the product species
``transitionState``    The string label of the transition state
``tunneling``          Method of estimating the quantum tunneling factor (optional)
====================== =========================================================

The following is an example of a typical reaction function::

    reaction(
        label = 'H + C2H4 <=> C2H5',
        reactants = ['H', 'C2H4'],
        products = ['C2H5'],
        transitionState = 'TS',
        tunneling='Eckart'        
    )

Note: the quantum tunneling factor method may be assigned either ``'Eckart'`` or ``'Wigner'``.

Thermodynamics Computations
===========================

Use a ``thermo()`` function to make CanTherm execute the thermodynamic
parameters computatiom for a species. Pass the string label of the species
you wish to compute the  thermodynamic parameters for and the type of
thermodynamics polynomial to generate (either ``'Wilhoit'`` or ''`NASA`'').
A table of relevant thermodynamic parameters will also be displayed in the
output file.

Below is a typical ``thermo()`` execution function::

    thermo('ethane', model='NASA')

Kinetics Computations
=====================

Use a ``kinetics()`` function to make CanTherm execute the high-pressure limit kinetic
parameters computation for a reaction. If desired, define a temperature range and number
of temperatures at which the high-pressure rate coefficient will be tabulated and saved to 
the outupt file. The 3-parameter modified Arrhenius coefficients will automatically be fit 
to the computed rate coefficients. The quantum tunneling factor will also be displayed.

Below is a typical ``kinetics()`` function::

    kinetics(    
    label = 'H + C2H4 <=> C2H5',
    Tmin = (400,'K'), Tmax = (1200,'K'), Tcount = 6,
    )

If specific temperatures are desired, you may specify a list
(``Tlist = ([400,500,700,900,1100,1200],'K')``) instead of Tmin, Tmax, and Tcount.

This is also acceptable::

    kinetics('H + C2H4 <=> C2H5')

Examples
========

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples`` directory.

Troubleshooting and FAQs
========================

1) The network that CanTherm generated and the resulting pdf file show abnormally large
absolute values. What's going on?

    This can happen if the number of atoms and atom types is not properly defined or consistent in your input file(s).

Cantherm User Checklist
========================

Using cantherm, or any rate theory package for that matter, requires careful consideration and management of a large amount of data, files, and input parameters. As a result, it is easy to make a mistake somewhere. This checklist was made to minimize such mistakes for users:

- Do correct paths exist for pointing to the files containing the electronic energies, molecular geometries and vibrational frequencies?

For calculations involving pressure dependence:

- Does the network pdf look reasonable? That is, are the relative energies what you expect based on the input?

For calculations using internal hindered rotors:

- Did you check to make sure the rotor has a reasonable potential (e.g., visually inspect the automatically generated rotor pdf files)?
- Within your input files, do all specified rotors point to the correct files?
- Do all of the atom label indices correspond to those in the file that is read by the logger (GaussianLog, QchemLog, etc.)?
- Why do the fourier fits look so much different than the results of the ab initio potential energy scan calculations? This is likely because the initial scan energy is not at a minimum. One solution is to simply shift the potential with respect to angle so that it starts at zero and, instead of having CanTherm read a Qchem or Gaussian output file, have CanTherm point to a 'ScanLog' file. Another problem can arise when the potential at 2*pi is also not [close] to zero.       
