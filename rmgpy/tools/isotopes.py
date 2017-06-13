#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains functionality for generating mechanisms with isotopes.
"""

import os
import os.path
import logging
import numpy as np
import itertools
from copy import copy, deepcopy
import pandas as pd
import shutil
import math

import rmgpy.constants as constants
from rmgpy.molecule import Molecule
from rmgpy.molecule.element import getElement
from rmgpy.tools.loader import loadRMGJob
from rmgpy.chemkin import ChemkinWriter
from rmgpy.rmg.input import getInput
from rmgpy.rmg.main import RMG, initializeLog
from rmgpy.rmg.model import Species
from rmgpy.species import Species as Species2
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.thermo.thermoengine import processThermoData
from rmgpy.data.thermo import findCp0andCpInf
from rmgpy.data.rmg import getDB

def initializeIsotopeModel(rmg, isotopes):
    """
    Initialize the RMG object by using the parameter species list
    as initial species instead of the species from the RMG input file.

    """
    # Read input file
    rmg.loadInput(rmg.inputFile)

    # Check input file 
    rmg.checkInput()

    # Load databases
    rmg.loadDatabase()

    logging.info("isotope: Adding the isotopomers into the RMG model")
    for spc in isotopes:
        spec, isNew = rmg.reactionModel.makeNewSpecies(spc)
        spec.thermo = spc.thermo        
        if isNew:
            rmg.reactionModel.addSpeciesToEdge(spec)
            rmg.initialSpecies.append(spec)
    logging.info("isotope: Adding standard species into the model")
    for spec in rmg.initialSpecies:
        spec.thermo = processThermoData(spec, spec.thermo)
        if not spec.reactive:
            rmg.reactionModel.enlarge(spec)
    for spec in rmg.initialSpecies:
        if spec.reactive:
            rmg.reactionModel.enlarge(spec)
    logging.info("isotope: Finalizing the species additions")
    rmg.initializeReactionThresholdAndReactFlags()
    rmg.reactionModel.initializeIndexSpeciesDict()


def generateIsotopeModel(outputDirectory, rmg0, isotopes):
    """
    Replace the core species of the rmg model with the parameter list
    of species.

    Generate all reactions between new list of core species.

    Returns created RMG object.
    """
    logging.debug("isotope: called generateIsotopeModel")
    rmg = RMG(inputFile=rmg0.inputFile, outputDirectory=outputDirectory)
    rmg.attach(ChemkinWriter(outputDirectory))

    logging.info("isotope: making the isotope model for with all species")
    initializeIsotopeModel(rmg, isotopes)

    logging.info("isotope: enlarging the isotope model")
    rmg.reactionModel.enlarge(reactEdge=True,
        unimolecularReact=rmg.unimolecularReact,
        bimolecularReact=rmg.bimolecularReact)

    logging.info("isotope: clustering reactions")
    clusters = cluster(rmg.reactionModel.core.reactions)
    logging.info('isotope: fixing the directions of every reaction to a standard')
    for isotopomerRxnList in clusters:
        ensureReactionDirection(isotopomerRxnList)

    consistent = True
    logging.info("isotope: checking symmetry is consistent amond isotopomers")
    for species_list in cluster(rmg.reactionModel.core.species):
        if not ensure_correct_symmetry(species_list):
            logging.info("isotopomers of {} with index {} may have wrong symmetry".format(species_list[0], species_list[0].index))
            consistent = False
    logging.info("isotope: checking that reaction degeneracy is consistent among isotopomers")
    for rxn_list in clusters:
        if not ensure_correct_degeneracies(rxn_list):
            logging.info("isotopomers of {} with index {} may have incorrect degeneracy.".format(rxn_list[0], rxn_list[0].index))
            consistent = False
    if not consistent:
        logging.warning("isotope: non-consistent degeneracy and/or symmetry was detected. This may lead to unrealistic deviations in enrichment. check log for more details")

    logging.info("isotope: saving files")
    rmg.saveEverything()

    rmg.finish()

    return rmg

def generateIsotopomers(spc, N=1):
    """
    Generate all isotopomers of the parameter species by adding max. N carbon isotopes to the
    atoms of the species.
    """

    mol = spc.molecule[0]
    isotope = getElement(6, 13)
    carbons = filter(lambda at: at.symbol == isotope.symbol, mol.atoms)

    mols = []
    addIsotope(0, N, mol, mols, isotope)

    spcs = []
    for isomol in mols:
        isotopomer = Species(molecule=[isomol], thermo=deepcopy(spc.thermo), transportData=spc.transportData, reactive=spc.reactive)
        isotopomer.generateResonanceIsomers(keepIsomorphic=True)
        spcs.append(isotopomer)

    # do not retain identical species:
    filtered = []
    while spcs:
        candidate = spcs.pop()
        unique = True
        for isotopomer in filtered:
            if isotopomer.isIsomorphic(candidate):
                unique = False
                break
        if unique: filtered.append(candidate)

    for isotopomer in filtered:
        correctEntropy(isotopomer, spc)

    return filtered

def addIsotope(i, N, mol, mols, element):
    """
    Iterate over the atoms of the molecule, and changes the element object
    of the atom by the provide parameter element object. Add the newly created
    isotopomer to the list of Molecule objects. For each created isotopomer,
    recursively call the method, until the maximum number of isotopes per molecule
    (N) is reached.

    """
    if i == N: return
    else:
        atoms = filter(lambda at: at.symbol == element.symbol, mol.atoms)
        for at in atoms:
            if at.element == element: continue
            else:
                isotopomer = mol.copy(deep=True)
                isotopomer.atoms[mol.atoms.index(at)].element = element
                mols.append(isotopomer)
                addIsotope(i+1, N, isotopomer, mols, element)

def cluster(objList):
    """
    Creates subcollections of isotopomers/reactions that 
    only differ in their isotopic labeling.

    This method works for either species or reactions.

    It is O(n^2) efficient
    """

    unclustered = copy(objList)

    # [[list of Species objs]]
    clusters = []

    while unclustered:
        candidate = unclustered.pop()
        for cluster in clusters:
            if compareIsotopomers(cluster[0],candidate):
                cluster.append(candidate)
                break
        else:
            clusters.append([candidate])

    return clusters

def removeIsotope(labeledObj, inplace = False):
    """
    Create a deep copy of the first molecule of the species object and replace
    non-normal Element objects (of special isotopes) by the 
    expected isotope.

    If the boolean `inplace` is True, the method remove the isotopic atoms of 
    the Species/Reaction
    inplace and returns a list of atom objects & element pairs for adding back
    to the oritinal object. This should significantly improve speed of this method.

    If successful, the non-inplace parts should be removed
    """

    if isinstance(labeledObj,Species2):
        if inplace:
            modifiedAtoms = []
            for mol in labeledObj.molecule:
                for atom in mol.atoms:
                    if atom.element.isotope != -1:
                        modifiedAtoms.append((atom,atom.element))
                        atom.element = getElement(atom.element.symbol)
            return modifiedAtoms
        else:
            stripped = labeledObj.copy(deep=True)
    
            for atom in stripped.molecule[0].atoms:
                if atom.element.isotope != -1:
                    atom.element = getElement(atom.element.symbol)

        # only do it for the first molecule, generate the other resonance isomers.
            stripped.molecule = [stripped.molecule[0]]
            stripped.generateResonanceIsomers()

        return stripped

    elif isinstance(labeledObj,Reaction):

        if inplace:

            atomList = []
            for reactant in  labeledObj.reactants:
                removed = removeIsotope(reactant,inplace)
                if removed:
                    atomList += removed
            for product in labeledObj.products:
                removed = removeIsotope(product,inplace)
                if removed:
                    atomList += removed

            return atomList
        else:
            strippedRxn = labeledObj.copy()

            strippedReactants = []
            for reactant in  strippedRxn.reactants:
                strippedReactants.append(removeIsotope(reactant,inplace))
            strippedRxn.reactants = strippedReactants

            strippedProducts = []
            for product in  strippedRxn.products:
                strippedProducts.append(removeIsotope(product,inplace))
            strippedRxn.products = strippedProducts

            return strippedRxn
    elif isinstance(labeledObj,Molecule):
        if inplace:
            modifiedAtoms = []
            for atom in labeledObj.atoms:
                if atom.element.isotope != -1:
                    modifiedAtoms.append((atom,atom.element))
                    atom.element = getElement(atom.element.symbol)
            return modifiedAtoms
        else:
            stripped = labeledObj.copy(deep=True)

            for atom in stripped.atoms:
                if atom.element.isotope != -1:
                    atom.element = getElement(atom.element.symbol)

            return stripped
    else:
        raise TypeError('Only Reaction, Species, and Molecule objects are supported')

def ensureReactionDirection(isotopomerRxns):
    """
    given a list of reactions with varying isotope labels but identical structure,
    obtained from the `cluster` method, this method remakes the kinetics so that 
    they all face the same direction.
    """

    # find isotopeless reaction as standard
    reference = isotopomerRxns[0]
    family = getDB('kinetics').families[reference.family]
    if family.ownReverse:
        for rxn in isotopomerRxns:
            if not compareIsotopomers(rxn, reference, eitherDirection=False):
                # the reaction is in the oposite direction
                logging.info('isotope: identified flipped reaction direction in reaction number {} of reaction {}. Altering the direction.'.format(rxn.index, str(rxn)))
                # reverse reactants and products
                rxn.reactants, rxn.products = rxn.products, rxn.reactants
                rxn.pairs = [(p,r) for r,p in rxn.pairs]

                # calculateDegeneracy
                rxnMols = TemplateReaction(reactants = [spec.molecule[0] for spec in rxn.reactants],
                                           products = [spec.molecule[0] for spec in rxn.products])
                forwardDegen = family.calculateDegeneracy(rxnMols)

                # set degeneracy to isotopeless reaction
                rxn.degeneracy = reference.degeneracy
                # make this reaction have kinetics of isotopeless reaction
                newKinetics = deepcopy(reference.kinetics)
                rxn.kinetics = newKinetics

                # set degeneracy to new reaction
                rxn.degeneracy = forwardDegen


def redoIsotope(atomList):
    """
    This takes a list of zipped atoms with their isotopes removed, from 
    and elements.
    """
    for atom, element in atomList:
        atom.element = element

def compareIsotopomers(obj1, obj2, eitherDirection = True):
    """
    This method takes two species or reaction objects and returns true if
    they only differ in isotopic labeling, and false if they have other
    differences. This also compares templates and families for TemplateReactions.

    The removeIsotope method can be slow, especially when comparing molecules
    and reactions. This was due to many copying of objects.

    This method avoid copying by storing the isotope and atom objects,
    removing them, doing the comparison, and rewriting them when
    finished the comparison.
    """

    atomlist = removeIsotope(obj1,inplace=True) + removeIsotope(obj2,inplace=True)
    if isinstance(obj1,Reaction):
        # make sure isotomorphic
        comparisonBool = obj1.isIsomorphic(obj2, eitherDirection)
        if comparisonBool and isinstance(obj1, TemplateReaction):
            # ensure families are the same
            comparisonBool = obj1.family == obj2.family
            if comparisonBool and not eitherDirection:
                # make sure templates are identical if in the same direction
                comparisonBool = frozenset(obj1.template) == frozenset(obj2.template)
    elif isinstance(obj1,Species2):
        comparisonBool = obj1.isIsomorphic(obj2)
    else:
        raise TypeError('Only Reaction and Speicies Objects are supported in compareIsotopomers')
    redoIsotope(atomlist)
    return comparisonBool

def generateRMGModel(inputFile, outputDirectory):
    """
    Generate the RMG-Py model NOT containing any non-normal isotopomers.

    Returns created RMG object.
    """
    initializeLog(logging.INFO, os.path.join(outputDirectory, 'RMG.log'))
    # generate mechanism:
    rmg = RMG(inputFile = os.path.abspath(inputFile),
            outputDirectory = os.path.abspath(outputDirectory)
        )
    rmg.execute()

    return rmg

def correctEntropy(isotopomer, isotopeless):
    """
    Correct the entropy of the isotopomer by the following correction for symmetry:

    S(corrected) = S(original) + R*ln(sigma(isotopeless)) - R*ln(sigma(isotopomer))

    This method also copies the Enthalpy, Cp and other thermo parameters from isotopeless
    """

    # calculate -R ln (sigma) in SI units (J/K/mol)
    Sisotopeless = - constants.R * math.log(isotopeless.getSymmetryNumber())
    Sisotopomer = - constants.R * math.log(isotopomer.getSymmetryNumber())

    # convert species thermo to ThermoData object:
    nasa = isotopomer.thermo

    # apply correction to entropy at 298K
    deltaS = Sisotopomer - Sisotopeless
    nasa = nasa.changeBaseEntropy(deltaS)

    # put the corrected thermo back as a species attribute:
    isotopomer.thermo = nasa

def isEnriched(obj):
    """
    Returns True if the species or reaction object has any enriched isotopes.
    """

    if isinstance(obj,Species):
        for atom in obj.molecule[0].atoms:
            if atom.element.isotope != -1 and not np.allclose(atom.element.mass, getElement(atom.element.symbol).mass):
                return True
        return False
    elif isinstance(obj,Reaction):
        enriched = []
        for spec in obj.reactants:
            enriched.append(isEnriched(spec))
        for spec in obj.products:
            enriched.append(isEnriched(spec))
        return any(enriched)
    else:
        raise TypeError('isEnriched only takes species and reaction objects. {} was sent'.format(str(type(obj))))


def ensure_correct_symmetry(isotopmoper_list, isotopic_element = 'C'):
    """
    given a list of isotopomers (species' objects) and the element that is labeled,
    returns True if the correct symmetry is detected. False if not detected

    This uses the observation that if there are less than 2^n isotopomors, where n is the number of
    isotopic elements, then there should be an 2^n - m isotopomers where 
    m is related to the amount of entropy increased by some isotopomers since they
    lost symmetry.
    """
    import math
    gas_constant = 8.314 # J/molK
    number_elements = 0
    for atom in isotopmoper_list[0].molecule[0].atoms:
        if atom.element.symbol ==isotopic_element:
            number_elements +=1

    minimum_entropy = min([spec.getEntropy(298) for spec in isotopmoper_list])

    count = 0.
    for spec in isotopmoper_list:
        entropy_diff = spec.getEntropy(298) - minimum_entropy
        count += math.exp(entropy_diff / gas_constant)
    return abs(count - 2**number_elements) < 0.01

def ensure_correct_degeneracies(reaction_isotopomer_list, print_data = False, r_tol_small_flux=1e-5, r_tol_deviation = 0.0001):
    """
    given a list of isotopomers (reaction objects), returns True if the correct 
    degeneracy values are detected. False if incorrect degeneracy values 
    exist. 

    This method assuumes an equilibrium distribution of compounds and finds 
    the fluxes created by the set of reactions (both forward and 
    reverse). It then checks that the fluxes of each compounds are proportional to their
    symmetry numbers (since lower symmetry numbers have higher entropy and should have a
    larger concentration.)

    inputs:

    reaction_isotopomer_list - a list of reactions that differ based on isotope placement
    print_data - output the table of fluxes obtained. useful for debugging incorrect values
    r_tol_small_flux - fraction of maximum flux to count as 'zero'
    r_tol_deviation - allowable numerical deviation in answers

    This method has a few datastructures it utilizes:

    product_structures - a list of tuples, (index, isotopomer_structure), containing reactants and products
    product_list - a pandas.DataFrame that sotres the fluxes and symmetry values
    """

    from rmgpy.kinetics.arrhenius import MultiArrhenius
    def store_flux_info(species, flux, product_list,  product_structures):
        """
        input:
            species - The desired species that you'd like to modify the flux value of
            flux - the amount to modify the flux of the species
            product_list - a list of current flux and symmetry values
            product_structures - a list of unlabled structures un the reaction class 

        modifies the product list by either adding on new row with the species' structure, flux
        symmetry ratio, etc. or adds the flux to the species' current row in `product_list`

        returns: modified product_list
        """
        # find the corresponding unlabeled structure
        for struc_index, product_structure in enumerate(product_structures):
            if compareIsotopomers(product_structure, species):
                structure_index = struc_index
                break 
        # store product flux and symmetry info
        for index in product_list.index:
            spec = product_list.at[index,'product']
            # if species already listed, add to its flux
            if product_list.at[index,'product_struc_index'] == structure_index \
                        and spec.isIsomorphic(species):
                product_list.at[index,'flux'] += flux
                return product_list
        else: # add product to list
            symmetry_ratio = product_structures[structure_index].getSymmetryNumber() / float(species.getSymmetryNumber())
            return product_list.append({'product': species, 
                                                'flux': flux,
                                                'product_struc_index': structure_index,
                                                'symmetry_ratio': symmetry_ratio},
                                               ignore_index=True)

    # copy list in case it is called upon
    reaction_isotopomer_list = copy(reaction_isotopomer_list)
    product_list = pd.DataFrame(columns=['product','flux','product_struc_index','symmetry_ratio'])
    product_structures = []
    unlabeled_rxn = None

    # check if first reaction is unlabeled. If not, make unlabeled reaction
    if not isEnriched(reaction_isotopomer_list[0]):
        unlabeled_rxn = reaction_isotopomer_list[0]
    else:
        unlabeled_rxn = reaction_isotopomer_list[0].copy()
        for mol in unlabeled_rxn.reactants:
            removeIsotope(mol, inplace = True)
        for mol in unlabeled_rxn.products:
            removeIsotope(mol, inplace = True)
    unlabeled_symmetry_reactants = np.prod([mol.getSymmetryNumber() for mol in unlabeled_rxn.reactants])
    unlabeled_symmetry_products = np.prod([mol.getSymmetryNumber() for mol in unlabeled_rxn.products])

    # prepare index of structures (product_structures)
    for struc in unlabeled_rxn.reactants + unlabeled_rxn.products:
        structure_index = len(product_structures)
        new_structure = struc.copy(deep=True)
        product_structures.append(new_structure)

    # go through each rxn cataloging fluxes to each product
    for rxn in reaction_isotopomer_list:
        # find characteristic flux for the forward direction
        reactant_conc = 1.
        reactant_conc /= np.prod([mol.getSymmetryNumber() for mol in rxn.reactants])
        reactant_conc *= unlabeled_symmetry_reactants
        if isinstance(rxn.kinetics, MultiArrhenius):
            rate = sum([arr.A.value_si for arr in rxn.kinetics.arrhenius])
        else:
            rate = rxn.kinetics.A.value_si
        product_flux = reactant_conc * rate

        # modify fluxes for forward direction
        for rxn_product in rxn.products:
            product_list = store_flux_info(rxn_product,product_flux, product_list, product_structures)
        for rxn_reactant in rxn.reactants:
            product_list = store_flux_info(rxn_reactant,-product_flux, product_list, product_structures)

        # now find characteristic flux of reverse direction. 
        if isinstance(rxn.kinetics, MultiArrhenius):
            reverse_A_factor = 0
            for arr in rxn.kinetics.arrhenius:
                reverse_A_factor += arr.A.value_si / rxn.getEquilibriumConstant(298)
        else:
            reverse_A_factor = rxn.kinetics.A.value_si / rxn.getEquilibriumConstant(298)

        # get reverse flux using product symmetries
        product_conc = 1.
        product_conc /= np.prod([mol.getSymmetryNumber() for mol in rxn.products])
        product_conc *= unlabeled_symmetry_products
        reactant_flux = product_conc * reverse_A_factor

        # modify reverse fluxes
        for rxn_product in rxn.products:
            product_list = store_flux_info(rxn_product,-reactant_flux, product_list, product_structures)
        for rxn_reactant in rxn.reactants:
            product_list = store_flux_info(rxn_reactant,reactant_flux, product_list, product_structures)

    if print_data:
        print product_list.sort_values(['product_struc_index','symmetry_ratio'])

    # now ensure the fluxes are correct or cancel out & throw error if not.
    pass_species = []
    for index, product_structure in enumerate(product_structures):
        products = product_list[product_list.product_struc_index == index]
        fluxes = products.flux / products.flux.sum()
        symmetries = products.symmetry_ratio / products.symmetry_ratio.sum()
        # the two decisison criteria
        accurate = np.allclose(fluxes,symmetries,rtol=r_tol_deviation)
        low_fluxes = all(products.flux.abs() < max(reactant_flux, product_flux)*r_tol_small_flux)
        pass_species.append(low_fluxes or accurate)
    return all(pass_species)

def run(inputFile, outputDir, original=None):
    """
    Accepts one input file with the RMG-Py model to generate. This file can contain 
    the optional parameter 'maximumIsotopicAtoms' to reduce isotope labeling.

    Firstly, generates the RMG model for the first input file. Takes the core species of that mechanism
    and generates all isotopomers of those core species. Next, generates all reactions between the
    generated pool of isotopomers, and writes it to file. 
    """
    logging.info("isotope: Starting the RMG isotope generation method 'run'")
    if not original:
        logging.info("isotope: original model not found, generating new one in directory `rmg`")
        logging.info("isotope: check `rmg/RMG.log` for the rest of the logging info.")

        outputdirRMG = os.path.join(outputDir, 'rmg')
        os.mkdir(outputdirRMG)

        rmg = generateRMGModel(inputFile, outputdirRMG)
    else:
        logging.info("isotope: original model being copied from previous RMG job in folder {}".format(original))
        outputdirRMG = original
        chemkinFile = os.path.join(outputdirRMG, 'chemkin', 'chem_annotated.inp')
        dictFile = os.path.join(outputdirRMG, 'chemkin', 'species_dictionary.txt')
        rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False, useChemkinNames=True)

    logging.info("isotope: generating isotope model")
    logging.info('Generating isotopomers for the core species in {}'.format(outputdirRMG))
    isotopes = []

    try:
        speciesConstraints = getInput('speciesConstraints')
    except Exception, e:
        logging.debug('Species constraints could not be found.')
        raise e

    try:
        maxIsotopes = speciesConstraints['maximumIsotopicAtoms']
    except KeyError, e:
        print ('Could not find the maxIsotopicAtoms constraint in the input file. Exiting...')
        raise e
    logging.info("isotope: adding all the new isotopomers")
    for spc in rmg.reactionModel.core.species:
        findCp0andCpInf(spc, spc.thermo)
        isotopes.extend(generateIsotopomers(spc, maxIsotopes))

    logging.info("isotope: adding original species to the model")
    # add the original unlabeled species:
    isotopes.extend(rmg.reactionModel.core.species)
    logging.info('Number of isotopomers: {}'.format(len(isotopes)))

    outputdirIso = os.path.join(outputDir, 'iso')
    os.mkdir(outputdirIso)

    logging.info('isotope: Generating RMG isotope model in {}'.format(outputdirIso))
    generateIsotopeModel(outputdirIso, rmg, isotopes)
