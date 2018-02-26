from mykrobe.probes import AlleleGenerator
from mykrobe.variant.schema import Variant
from mykrobe.variant.schema import VariantSet
from mykrobe.variant.schema import Reference
from mykrobe.variant.schema import ReferenceSet
from nose.tools import assert_raises
from mongoengine import connect
DB = connect('atlas-test')


class TestINDELAlleleGenerator():

    def setup(self):
        DB.drop_database('atlas-test')

        self.pg = AlleleGenerator(
            reference_filepath="atlasvar/data/BX571856.1.fasta")
        self.pg2 = AlleleGenerator(
            reference_filepath="atlasvar/data/NC_000962.2.fasta")
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file",
            reference_set=self.reference_set)
        self.variant_sets = [self.variant_set]
        self.reference = Reference().create_and_save(
            name="ref",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])

    def test_simple_deletion1(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AA",
            start=31,
            alternate_bases=["A"])
        assert v.is_indel
        assert v.is_deletion
        panel = self.pg.create(v)
        assert panel.ref == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert self.pg._calculate_length_delta_from_indels(v, []) == 1
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTGATC"]

    def test_simple_deletion2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AT",
            start=32,
            alternate_bases=["A"])
        panel = self.pg.create(v)
        assert panel.ref == "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == [
            "GATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGATCC"]

    def test_simple_deletion3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AT",
            start=2902618,
            alternate_bases=["T"])
        panel = self.pg.create(v)
        assert panel.ref == "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert panel.alts == [
            "ATAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTT"]

    def test_simple_deletion4(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"])
        panel = self.pg.create(v)
        assert panel.ref == "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"]

    def test_simple_insertion1(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=1,
            alternate_bases=["TTTC"])
        panel = self.pg.create(v)
        assert v.is_indel
        assert v.is_insertion
        assert panel.ref == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == [
            "TTTCGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=1,
            alternate_bases=["CTTT"])
        panel = self.pg.create(v)
        assert panel.ref == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == [
            "CTTTGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["ATTT"])
        panel = self.pg.create(v)
        assert panel.ref == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion4(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["AGGGG"])
        panel = self.pg.create(v)
        assert panel.ref == "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == [
            "GATTAAAGATAGAAATACACGATGCGAGCAAGGGGTCAAATTTCATAACATCACCATGAGTTTGATC"]

    def test_simple_insertion5(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=2902618,
            alternate_bases=["ATGC"])
        panel = self.pg.create(v)
        assert panel.ref == "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert panel.alts == [
            "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTATGCT"]

    def test_double_insertion(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=4021408,
            alternate_bases=["ACGCTGGCGGGCG"])
        v1 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AGA",
            start=4021406,
            alternate_bases=["CGG"])
        context = [v1]
        assert self.pg2._remove_overlapping_contexts(v, [v1]) == []
        panel = self.pg2.create(v, context=context)
        assert panel.ref == "ATCTAGCCGCAAGGGCGCGAGCAGACGCAGAATCGCATGATTTGAGCTCAAATCATGCGATTC"
        assert panel.alts == [
            "ATCTAGCCGCAAGGGCGCGAGCAGACGCAGACGCTGGCGGGCGATCGCATGATTTGAGCTCAAATCATGCGATTC"]

    def test_large_insertion(self):
        v = Variant.create(variant_sets=self.variant_sets, reference=self.reference, reference_bases="CCGCCGGCCCCGCCGTTT", start=1636155, alternate_bases=[
                           "CTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCG"])
        panel = self.pg2.create(v, context=[])
        assert panel.ref == "GAGTCGCCGAGGACGCCGGCGCCGCCATTGTCGCCAAATACCGTGAGACCTAGCAGGGTGCCGGCGCCGCCCTTGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCATCGCCGATGATGTTTTCCCCGCCCTTGCCGCCAGCCCCAGCGTTCCCG"
        assert panel.alts == [
            "GGTTGGATCGCCACCGGCGCCACCGGCGCCGCCCGCGCCACCAGCACCGCCGCTGCCATCTGGGTCCGTCGAGTCGCCGAGGACGCCGGCGCCGCCATTGTCGCCAAATACCGTGAGACCTAGCAGGGTGCCGGCGCCGCCCTTGCTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCGCCATCGCCGATGATGTTTTCCCCGCCCTTGCCGCCAGCCCCAGCGTTCCCGCCGGCTCCGCCACTGGCGCCGGTGCCGCCGGGTGCAACGGCGTTGGCGCCGTTACCGCCGTTGCCGCCTTT"]
