# -*- coding: utf-8 -*-

"""Tests for conversion utilities for BEL."""

import unittest
from typing import Mapping

from pybel import BELGraph
from pybel.constants import NEGATIVE_CORRELATION, POSITIVE_CORRELATION, REGULATES
from pybel.dsl import Abundance, BaseEntity, Gene, Pathology, Protein, Reaction, activity, pmod, secretion
from pybel.testing.utils import n
from pybel_obo import convert

# the n() function makes a dummy name
abundance_1 = Abundance(namespace=n(), identifier=n(), name=n())
pathology_1 = Pathology(namespace=n(), identifier=n(), name=n())
protein_1 = Protein(namespace=n(), identifier=n(), name=n())
protein_2 = Protein(namespace=n(), identifier=n(), name=n())
rna_1 = protein_1.get_rna()
rna_2 = protein_2.get_rna()
gene_1 = rna_1.get_gene()
gene_2 = Gene(namespace=n(), identifier=n(), name=n())


def relationship(ro_id: str, label: str) -> Mapping[str, str]:
    """Build a relationship dictionary."""
    return dict(ro_id=ro_id, label=label)


class TestConversion(unittest.TestCase):
    """Test several conversions."""

    def setUp(self) -> None:
        """Set up each test with an instance of a BEL graph."""
        self.graph = BELGraph()

    def convert(self, source: BaseEntity, target: BaseEntity, key: str) -> Mapping[str, str]:
        """Wrap conversion on this test case's graph."""
        return convert(self.graph, source, target, key)

    def test_regulates_transport(self):
        """Test generating a `regulates transport of` (RO:0002011) relationship."""
        r = relationship('RO:0002011', 'regulates transport of')

        key = self.graph.add_increases(protein_1, protein_2, n(), n(), object_modifier=secretion())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

        key = self.graph.add_directly_increases(protein_1, protein_2, n(), n(), object_modifier=secretion())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

        key = self.graph.add_decreases(protein_1, protein_2, n(), n(), object_modifier=secretion())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

        key = self.graph.add_directly_decreases(protein_1, protein_2, n(), n(), object_modifier=secretion())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

        key = self.graph.add_qualified_edge(protein_1, protein_2, REGULATES, n(), n(), object_modifier=secretion())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

    def test_activity_directly_regulates_activity_of(self):
        """Test generating an `activity directly regulates activity of` (RO:0002448) relationship."""
        r = relationship('RO:0002448', 'activity directly regulates activity of')

        key = self.graph.add_qualified_edge(protein_1, protein_2, REGULATES, n(), n(),
                                            subject_modifier=activity(),
                                            object_modifier=activity())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

    def test_activity_directly_negatively_regulates_activity_of(self):
        """Test generating an `activity directly negatively regulates activity of` (RO:0002449) relationship."""
        r = relationship('RO:0002449', 'activity directly negatively regulates activity of')

        key = self.graph.add_directly_decreases(protein_1, protein_2, n(), n(),
                                                subject_modifier=activity(),
                                                object_modifier=activity())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

    def test_represses_expression_of(self):
        """Test generating a `represses expression of` (RO:0003002) relationship."""
        r = relationship('RO:0003002', 'represses expression of')

        key = self.graph.add_directly_decreases(protein_1, rna_2, n(), n(),
                                                subject_modifier=activity())
        self.assertEqual(r, self.convert(protein_1, rna_2, key))

        # activity of the subject isn't strictly necessary
        key = self.graph.add_directly_decreases(protein_1, rna_2, n(), n())
        self.assertEqual(r, self.convert(protein_1, rna_2, key))

    def test_activity_directly_positively_regulates_activity_of(self):
        """Test generating an `activity directly positively regulates activity of` (RO:0002449) relationship."""
        r = relationship('RO:0002450', 'activity directly positively regulates activity of')

        key = self.graph.add_directly_increases(protein_1, protein_2, n(), n(),
                                                subject_modifier=activity(),
                                                object_modifier=activity())
        self.assertEqual(r, self.convert(protein_1, protein_2, key))

    def test_increases_expression_of(self):
        """Test generating an `increases expression of` (RO:0003002) relationship."""
        r = relationship('RO:0003003', 'increases expression of')

        key = self.graph.add_directly_increases(protein_1, rna_2, n(), n(),
                                                subject_modifier=activity())
        self.assertEqual(r, self.convert(protein_1, rna_2, key))

        # activity of the subject isn't strictly necessary
        key = self.graph.add_directly_increases(protein_1, rna_2, n(), n())
        self.assertEqual(r, self.convert(protein_1, rna_2, key))

    def test_is_substance_that_treats(self):
        """Test generating an `is substance that treats` (RO:0002606) relationship."""
        r = relationship('RO:0002606', 'is substance that treats')

        key = self.graph.add_decreases(abundance_1, pathology_1, n(), n())
        self.assertEqual(r, self.convert(abundance_1, pathology_1, key))

    def test_causes_condition(self):
        """Test generating a `causes condition` (RO:0003303) relationship."""
        r = relationship('RO:0003303', 'causes condition')

        key = self.graph.add_increases(abundance_1, pathology_1, n(), n())
        self.assertEqual(r, self.convert(abundance_1, pathology_1, key))

        key = self.graph.add_increases(gene_2, pathology_1, n(), n())
        self.assertEqual(r, self.convert(gene_2, pathology_1, key))

        key = self.graph.add_directly_increases(gene_2, pathology_1, n(), n())
        self.assertEqual(r, self.convert(gene_2, pathology_1, key))

    def test_over_expressed_in(self):
        """Test generating an `over-expressed in` (RO:0002245) relationship."""
        r = relationship('RO:0002245', 'over-expressed in')
        key = self.graph.add_qualified_edge(protein_1, pathology_1, POSITIVE_CORRELATION, n(), n())
        self.assertEqual(r, self.convert(protein_1, pathology_1, key))

    def test_under_expressed_in(self):
        """Test generating an `under-expressed in` (RO:0002246) relationship."""
        r = relationship('RO:0002246', 'under-expressed in')
        key = self.graph.add_qualified_edge(protein_1, pathology_1, NEGATIVE_CORRELATION, n(), n())
        self.assertEqual(r, self.convert(protein_1, pathology_1, key))

    def test_transcribed_from(self):
        """Test generating a `transcribed to` (RO:0002511) relationship."""
        r = relationship('RO:0002511', 'transcribed to')
        key = self.graph.add_translation(rna_1, protein_1)
        self.assertEqual(r, self.convert(rna_1, protein_1, key))

    def test_ribosomally_translates_to(self):
        """Test generating a `ribosomally translates to` (RO:0002513) relationship."""
        r = relationship('RO:0002513', 'ribosomally translates to')
        key = self.graph.add_transcription(gene_1, rna_1)
        self.assertEqual(r, self.convert(gene_1, rna_1, key))

        # TODO should we directly infer gene product of (RO:0002204)?

    def test_phosphorylates(self):
        """Test generating a `phosphorylates` (RO:0002447) relationship."""
        r = relationship('RO:0002447', 'phosphorylates')
        v_phosphorylated = protein_2.with_variants(pmod('Ph'))
        v_phosphorylation_reaction = Reaction(protein_2, v_phosphorylated)

        # Check simple BEL-style relations
        key = self.graph.add_increases(protein_1, v_phosphorylated, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_phosphorylated, key))

        key = self.graph.add_directly_increases(protein_1, v_phosphorylated, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_phosphorylated, key))

        # Check mechanistic BioPAX-style relations
        key = self.graph.add_increases(protein_1, v_phosphorylation_reaction, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_phosphorylation_reaction, key))

        key = self.graph.add_directly_increases(protein_1, v_phosphorylation_reaction, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_phosphorylation_reaction, key))

    def test_ubiquitinates(self):
        """Test generating a `ubiquitinates` (RO:0002480) relationship."""
        r = relationship('RO:0002480', 'ubiquitinates')
        v_ubiquitinated = protein_2.with_variants(pmod('Ub'))
        v_ubiquitination_reaction = Reaction(protein_2, v_ubiquitinated)

        # Check simple BEL-style relations
        key = self.graph.add_increases(protein_1, v_ubiquitinated, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_ubiquitinated, key))

        key = self.graph.add_directly_increases(protein_1, v_ubiquitinated, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_ubiquitinated, key))

        # Check mechanistic BioPAX-style relations
        key = self.graph.add_increases(protein_1, v_ubiquitination_reaction, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_ubiquitination_reaction, key))

        key = self.graph.add_directly_increases(protein_1, v_ubiquitination_reaction, n(), n())
        self.assertEqual(r, self.convert(protein_1, v_ubiquitination_reaction, key))

        # TODO what about other PTMs?
