import unittest
import os
import shutil
import numpy as np

from rmgpy.tools.loader import loadRMGJob
import rmgpy
from rmgpy.species import Species
from rmgpy.tools.isotopes import *
from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics.arrhenius import Arrhenius

from rmgpy.data.rmg import RMGDatabase

def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.load(
        path=os.path.join(settings['test_data.directory'], 'testing_database'),
        kineticsFamilies=[
            'H_Abstraction',
        ],
        testing=True,
        depository=False,
        solvation=False,
    )
    database.loadForbiddenStructures()

    # Prepare the database by loading training reactions and averaging the rate rules
    for family in database.kinetics.families.values():
        family.addKineticsRulesFromTrainingSet(thermoDatabase=database.thermo)
        family.fillKineticsRulesByAveragingUp(verbose=True)


def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg
    rmg.database = None


class IsotopesTest(unittest.TestCase):

    def testClusterWithSpecies(self):
        """
        Test that isotope partitioning algorithm work with Reaction Objects.
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        meth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
        """
    )

        spcList = [eth, ethi]

        clusters = cluster(spcList)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        spcList = [meth, ethi]

        clusters = cluster(spcList)
        self.assertEquals(len(clusters), 2)
        self.assertEquals(len(clusters[0]), 1)

    def testClusterWithReactions(self):
        """
        Test that isotope partitioning algorithm works with Reaction objects
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        meth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
        """
    )

        rxn0 = Reaction(reactants=[ethi,ethi], products=[ethi,eth])
        rxn1 = Reaction(reactants=[eth,ethi], products=[eth,eth])
        rxn2 = Reaction(reactants=[ethi,meth], products=[meth,ethi])
        rxn3 = Reaction(reactants=[eth,meth], products=[eth,meth])
        rxn4 = Reaction(reactants=[meth], products=[meth])
        rxn5 = Reaction(reactants=[ethi], products=[eth])
        
        
        sameClusterList = [rxn0,rxn1]

        clusters = cluster(sameClusterList)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)
        
        sameClusterList = [rxn2,rxn3]

        clusters = cluster(sameClusterList)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        multiClusterList = [rxn0, rxn1, rxn2, rxn3, rxn4, rxn5]

        clusters = cluster(multiClusterList)
        self.assertEquals(len(clusters), 4)
        self.assertEquals(len(clusters[0]), 1)

        
    def testRemoveIsotopeForReactions(self):
        """
        Test that remove isotope algorithm works with Reaction objects.
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )
        unlabeledRxn = Reaction(reactants=[eth], products = [eth])
        labeledRxn = Reaction(reactants=[ethi], products = [ethi])
        stripped = removeIsotope(labeledRxn)
        
        self.assertTrue(unlabeledRxn.isIsomorphic(stripped))

    def testRemoveIsotopeForSpecies(self):
        """
        Test that remove isotope algorithm works with Species.
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        stripped = removeIsotope(ethi)

        self.assertTrue(eth.isIsomorphic(stripped))

    def testInplaceRemoveIsotopeForReactions(self):
        """
        Test that removeIsotope and redoIsotope works with reactions
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )
        unlabeledRxn = Reaction(reactants=[eth], products = [eth])
        labeledRxn = Reaction(reactants=[ethi], products = [ethi])
        storedLabeledRxn = labeledRxn.copy()
        modifiedAtoms = removeIsotope(labeledRxn, inplace=True)

        self.assertTrue(unlabeledRxn.isIsomorphic(labeledRxn))

        redoIsotope(modifiedAtoms)

        self.assertTrue(storedLabeledRxn.isIsomorphic(labeledRxn))

    def testCompareIsotopomersWorksOnSpecies(self):
        """
        Test that compareIsotomers works on species objects
        """
        ethii = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """)
        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """)
        self.assertTrue(compareIsotopomers(ethii,ethi))

    def testCompareIsotopomersDoesNotAlterSpecies(self):
        """
        Test that compareIsotomers works on species objects
        """
        ethii = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """)
        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """)
        compareIsotopomers(ethii,ethi)
        
        # ensure species still have labels
        for atom in ethii.molecule[0].atoms:
            if atom.element.symbol == 'C':
                self.assertEqual(atom.element.isotope, 13, 'compareIsotopomer removed the isotope of a species.')

    def testCompareIsotopomersFailsOnSpecies(self):
        """
        Test that compareIsotomers fails on species objects
        """
        ethane = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """)
        ethenei = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 i13 {1,D} {6,S} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
        """)
        self.assertFalse(compareIsotopomers(ethane,ethenei))

    def testCompareIsotopomersWorksOnReactions(self):
        """
        Test that compareIsotomers works on different reaction objects
        """
        h = Species().fromAdjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """)
        h2 = Species().fromAdjacencyList(
        """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """)
        propanei = Species().fromAdjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 i13 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """)
        propane = Species().fromAdjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """)
        npropyli = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 i13 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """)
        npropyl = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """)

        reaction2 = TemplateReaction(reactants = [propanei,h],
                            products = [npropyli,h2],
                            family = 'H_Abstraction')
        reaction3 = TemplateReaction(reactants = [propane,h],
                            products = [h2,npropyl],
                            family = 'H_Abstraction')
        self.assertTrue(compareIsotopomers(reaction2,reaction3))
        
    def testCompareIsotopomersFailsOnReactions(self):
        """
        Test that compareIsotomers fails on different reaction objects
        """
        h = Species().fromAdjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """)
        h2 = Species().fromAdjacencyList(
        """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """)
        propanei = Species().fromAdjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 i13 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """)
        propane = Species().fromAdjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """)
        npropyli = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 i13 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """)
        npropyl = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """)

        reaction2 = TemplateReaction(reactants = [propanei,h],
                            products = [npropyli,h2],
                            family = 'H_Abstraction')
        
        magicReaction = TemplateReaction(reactants = [propane,h],
                            products = [propanei,h],
                            family = 'H_Abstraction')
        self.assertFalse(compareIsotopomers(reaction2,magicReaction))

    def testGenerateIsotopomers(self):
        """
        Test that the generation of isotopomers with N isotopes works.
        """
        from rmgpy.thermo.nasa import NASAPolynomial, NASA

        spc = Species().fromSMILES('CC')

        polynomial = NASAPolynomial(coeffs=[1.,1.,1.,1.,1.,1.,1.],
                         Tmin=(200,'K'),Tmax=(1600,'K'),E0=(1.,'kJ/mol'),
                         comment='made up thermo')
        
        spc.thermo = NASA(polynomials=[polynomial],Tmin=(200,'K'),
                        Tmax=(1600,'K'),E0=(1.,'kJ/mol'),
                        comment='made up thermo')

        spcs = generateIsotopomers(spc, 0)
        self.assertEquals(len(spcs), 0)

        spcs = generateIsotopomers(spc)
        self.assertEquals(len(spcs), 1)

        spcs = generateIsotopomers(spc, 2)
        self.assertEquals(len(spcs), 2)

        spcs = generateIsotopomers(spc, 3)
        self.assertEquals(len(spcs), 2)

    def testIsEnriched(self):
        """
        ensures the Enriched method functions
        """
        npropyl = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """)
        
        npropyli = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 i13 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """)
        
        self.assertTrue(isEnriched(npropyli))
        
        self.assertFalse(isEnriched(npropyl))
        
        enrichedReaction = TemplateReaction(reactants = [npropyl],
                            products = [npropyli],
                            family = 'H_Abstraction')
        self.assertTrue(isEnriched(enrichedReaction))
               
        bareReaction = TemplateReaction(reactants = [npropyl],
                            products = [npropyl],
                            family = 'H_Abstraction')
                            
        self.assertFalse(isEnriched(bareReaction))
