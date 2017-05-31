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
    differences.

    The removeIsotope method can be slow, especially when comparing molecules
    and reactions. This was due to many copying of objects.

    This method avoid copying by storing the isotope and atom objects,
    removing them, doing the comparison, and rewriting them when
    finished the comparison.
    """

    atomlist = removeIsotope(obj1,inplace=True) + removeIsotope(obj2,inplace=True)
    if isinstance(obj1,Reaction):
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
