from mykrobe.probes import AlleleGenerator
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantSet
from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet
from mongoengine import connect
DB = connect('mykrobe-test')
from base import assert_no_overlapping_snps
from mykrobe.utils import seq_to_kmers


class TestINDELAlleleGenerator():

    def setup(self):
        DB.drop_database('mykrobe-test')

        self.pg = AlleleGenerator(
            reference_filepath="src/mykrobe/data/BX571856.1.fasta")
        self.pg2 = AlleleGenerator(
            reference_filepath="src/mykrobe/data/NC_000962.3.fasta")
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
        assert_no_overlapping_snps(panel)                     
        assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
        assert self.pg._calculate_length_delta_from_indels(v, []) == 1
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTGA"]

    def test_simple_deletion2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AT",
            start=32,
            alternate_bases=["A"])
        panel = self.pg.create(v)
        assert_no_overlapping_snps(panel)                     
        assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
        assert panel.alts == [
            "GATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_deletion3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AT",
            start=2902618,
            alternate_bases=["T"])
        panel = self.pg.create(v)
        assert_no_overlapping_snps(panel)                     
        assert "TTTATACTACTGCTCAATTTTTTTACTTTTATNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" in panel.refs
        assert panel.alts == [
            "TTTATACTACTGCTCAATTTTTTTACTTTTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"]

    def test_simple_deletion4(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"])
        panel = self.pg.create(v)
        assert_no_overlapping_snps(panel)                     
        assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion1(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=1,
            alternate_bases=["TTTC"])
        panel = self.pg.create(v)
#        assert_no_overlapping_snps(panel)### Skip this test for vars in first k bases of ref                    
        assert v.is_indel
        assert v.is_insertion
        assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
        assert panel.alts == [
            "TTTCGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"]

    def test_simple_insertion2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=1,
            alternate_bases=["CTTT"])
        panel = self.pg.create(v)
#        assert_no_overlapping_snps(panel)### Skip this test for vars in first k bases of ref                    
        assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
        assert panel.alts == [
            "CTTTGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"]

    def test_simple_insertion3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["ATTT"])
        panel = self.pg.create(v)
        assert_no_overlapping_snps(panel)                     
        assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTG"]

    def test_simple_insertion4(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["AGGGG"])
        panel = self.pg.create(v)
        assert_no_overlapping_snps(panel)                     
        assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
        assert panel.alts == [
            "GATTAAAGATAGAAATACACGATGCGAGCAAGGGGTCAAATTTCATAACATCACCATGAGTTTGA"]

    def test_simple_insertion5(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=2902618,
            alternate_bases=["ATGC"])
        panel = self.pg.create(v)
        assert_no_overlapping_snps(panel)                     
        assert "TTTATACTACTGCTCAATTTTTTTACTTTTATNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" in panel.refs
        assert panel.alts == [
            "TTTATACTACTGCTCAATTTTTTTACTTTTATGCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"]

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
        assert_no_overlapping_snps(panel)                     
        assert "ATCTAGCCGCAAGGGCGCGAGCAGACGCAGAATCGCATGATTTGAGCTCAAATCATGCGAT" in panel.refs
        assert panel.alts == [
            "ATCTAGCCGCAAGGGCGCGAGCAGACGCAGACGCTGGCGGGCGATCGCATGATTTGAGCTCAAATCATGCGAT"]

    def test_double_indel_fail(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="CCA",
            start=2288851,
            alternate_bases=["A"])
        v1 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=2288850,
            alternate_bases=["ACC"])
        context = [v1]
        panel = self.pg2.create(v, context=context)
        assert_no_overlapping_snps(panel)                     
        assert "GGCGCACACAATGATCGGTGGCAATACCGACCACATCGACCTCATCGACGCCGCGTTGCCG" in panel.refs
        assert "GGCGCACACAATGATCGGTGGCAATACCGACCACATCGACCTCATCGACGCCGCGTTGCCG" not in panel.alts 

    def test_large_insertion(self):
        v = Variant.create(variant_sets=self.variant_sets, reference=self.reference, reference_bases="CCGCCGGCCCCGCCGTTT", start=1636155, alternate_bases=[
                           "CTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCG"])
        panel = self.pg2.create(v, context=[])
        assert_no_overlapping_snps(panel)                     
        assert "AGACCTAGCAGGGTGCCGGCGCCGCCCTTGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCAT" in panel.refs
        assert panel.alts == [
            "GGTTGGATCGCCACCGGCGCCACCGGCGCCGCCCGCGCCACCAGCACCGCCGCTGCCATCTGGGTCCGTCGAGTCGCCGAGGACGCCGGCGCCGCCATTGTCGCCAAATACCGTGAGACCTAGCAGGGTGCCGGCGCCGCCCTTGCTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCGCCATCGCCGATGATGTTTTCCCCGCCCTTGCCGCCAGCCCCAGCGTTCCCGCCGGCTCCGCCACTGGCGCCGGTGCCGCCGGGTGCAACGGCGTTGGCGCCGTTACCGCCGTTGCCGCCT"]
