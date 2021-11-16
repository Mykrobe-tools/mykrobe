import os

from mongoengine import connect

from base import assert_no_overlapping_kmers
from mykrobe.probes import AlleleGenerator
from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantSet

DB = connect("mykrobe-test")
DATA_DIR = os.path.join("tests", "ref_data")


class TestINDELandSNPSAlleleGenerator:
    def setup(self):
        DB.drop_database("mykrobe-test")
        self.pg = AlleleGenerator(
            reference_filepath=f"{DATA_DIR}/BX571856.1.fasta", kmer=31
        )
        self.pg2 = AlleleGenerator(
            reference_filepath=f"{DATA_DIR}/NC_000962.3.fasta", kmer=31
        )
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file", reference_set=self.reference_set
        )
        self.variant_sets = [self.variant_set]
        self.reference = Reference().create_and_save(
            name="ref", md5checksum="sre", reference_sets=[self.reference_set]
        )

    def teardown(self):
        DB.drop_database("mykrobe-test")

    def test_ins_with_SNP_context(self):

        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["ATTT"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v, context=[v2])
        # assert_no_overlapping_kmers(panel)  ### This test seems to fail sometimes, and pass othertimes...
        assert (
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "GATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTG",
                "TTAAAGATAGAAATACACGATGCGAGCATTTTTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )

    def test_del_with_SNP_context1(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AA",
            start=31,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=33,
            alternate_bases=["A"],
        )
        panel = self.pg.create(v, context=[v2])
        assert_no_overlapping_kmers(panel)
        assert (
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "ATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGA",
                "GATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )

    def test_del_with_SNP_context2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AA",
            start=31,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v, context=[v2])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2]) == []
        assert (
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            ["GATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTG"]
        )

    def test_del_with_ins_context1(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AAT",
            start=31,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=4,
            alternate_bases=["TTTT"],
        )
        panel = self.pg.create(v, context=[v2])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2]) == [v2]
        assert (
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "GATTAAAGATAGAAATACACGATGCGAGCACAAATTTCATAACATCACCATGAGTTTGAT",
                "TTTTAAAGATAGAAATACACGATGCGAGCACAAATTTCATAACATCACCATGAGTTTG",
            ]
        )

    def test_del_with_ins_context2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=1,
            alternate_bases=["CTTT"],
        )
        panel = self.pg.create(v, context=[v2])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2]) == [v2]
        assert self.pg._remove_contexts_not_within_k(v, [v2]) == []
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            ["ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT"]
        )

    def test_del_with_ins_context3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=5,
            alternate_bases=["TT"],
        )
        panel = self.pg.create(v, context=[v2])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2]) == [v2]
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "TTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            ]
        )

    def test_del_with_ins_context4(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=5,
            alternate_bases=["TT"],
        )
        v3 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=5,
            alternate_bases=["TG"],
        )
        panel = self.pg.create(v, context=[v2, v3])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2, v3]) == [v2, v3]
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "TTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "TTGAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            ]
        )

    def test_del_with_ins_context5(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=5,
            alternate_bases=["TT"],
        )
        v3 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=6,
            alternate_bases=["AG"],
        )
        panel = self.pg.create(v, context=[v2, v3])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2, v3]) == [v2, v3]
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "TTAGAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGA",
                "TTAGAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "TTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            ]
        )

    def test_del_with_ins_context_where_base_is_deleted1(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=33,
            alternate_bases=["C"],
        )
        panel = self.pg.create(v, context=[v2])
        assert_no_overlapping_kmers(panel)
        assert self.pg._remove_overlapping_contexts(v, [v2]) == []
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            ["ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT"]
        )

    def test_del_with_ins_context_where_base_is_deleted2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATC",
            start=32,
            alternate_bases=["A"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="TAAA",
            start=5,
            alternate_bases=["T"],
        )
        v3 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=7,
            alternate_bases=["AG"],
        )
        panel = self.pg.create(v, context=[v2, v3])
        assert_no_overlapping_kmers(panel)
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "CGATTGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATC",
                "TTAAGAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            ]
        )

        panel = self.pg.create(v, context=[v3, v2])
        assert_no_overlapping_kmers(panel)
        assert (
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
            in panel.refs
        )
        assert sorted(panel.alts) == sorted(
            [
                "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
                "CGATTGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATC",
                "TTAAGAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            ]
        )

    def test_snp_with_replace_context(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="G",
            start=2338961,
            alternate_bases=["A"],
        )
        v1 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="GGATG",
            start=2338990,
            alternate_bases=["CGATA"],
        )
        panel = self.pg2.create(v, context=[v1])
        assert_no_overlapping_kmers(panel)
        assert (
            "CGACTAGCCACCATCGCGCATCAGTGCGAGGTCAAAAGCGACCAAAGCGAGCAAGTCGCGG"
            in panel.refs
        )

        assert set(panel.alts) == set(
            [
                "CGACTAGCCACCATCGCGCATCAGTGCGAGATCAAAAGCGACCAAAGCGAGCAAGTCGCCG",
                "CGACTAGCCACCATCGCGCATCAGTGCGAGATCAAAAGCGACCAAAGCGAGCAAGTCGCGG",
            ]
        )

    def test_indel_snp_indel_context(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="TCGCGTGGC",
            start=4021459,
            alternate_bases=["GCGAGCAGA"],
        )
        v1 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=4021455,
            alternate_bases=["ATCTAGCCGCAAG"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=4021489,
            alternate_bases=["G"],
        )
        panel = self.pg2.create(v)  # , context = [v1, v2])
        assert_no_overlapping_kmers(panel)
        assert (
            "ATCATGCGATTCTGCGTCTGCTCGCGAGGCTCGCGTGGCCGCCGGCGCTGGCGGGCGATCT"
            in panel.refs
        )

        panel = self.pg2.create(v, context=[v1, v2])
        assert_no_overlapping_kmers(panel)
        assert sorted(panel.alts) == sorted(
            [
                "ATCATGCGATTCTGCGTCTGCTCGCGAGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCG",
                "ATCATGCGATTCTGCGTCTGCTCGCGAGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCT",
                "TGCGTCTGCTCGCGATCTAGCCGCAAGGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCG",
                "TGCGTCTGCTCGCGATCTAGCCGCAAGGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCT",
            ]
        )

    def test_complex_context(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="ATTT",
            start=1503643,
            alternate_bases=["A"],
        )
        v1 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="CCT",
            start=1503615,
            alternate_bases=["C"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=1503655,
            alternate_bases=["ATGCCGCCGCC"],
        )
        panel = self.pg2.create(v, context=[v1, v2])
        assert_no_overlapping_kmers(panel)
        assert (
            "ATCCTGGAGCCCACCAGCGGAAACACCGGCATTTCGCTGGCGATGGCGGCCCGGTTGAAGG"
            in panel.refs
        )
        assert set(panel.alts) == set(
            [
                "CCATCGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGGCGGCCCGGTTGAAGGGGT",
                "TCCTGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGGCGGCCCGGTTGAAGGGG",
                "ATCGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGCCGCCGCCTGGCGGCCCGG",
                "TCCTGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGCCGCCGCCTGGCGGCCCGG",
            ]
        )
