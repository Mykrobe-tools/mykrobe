import os

import pytest
from mongoengine import connect

from base import assert_no_overlapping_kmers
from mykrobe.probes import AlleleGenerator
from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantSet

DB = connect("mykrobe-test")
DATA_DIR = os.path.join("tests", "ref_data")


class TestSNPAlleleGenerator:
    def setup(self):
        DB.drop_database("mykrobe-test")
        self.pg = AlleleGenerator(
            reference_filepath=f"{DATA_DIR}/BX571856.1.fasta", kmer=31
        )
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file", reference_set=self.reference_set
        )
        self.variant_sets = [self.variant_set]
        self.reference = Reference().create_and_save(
            name="ref", md5checksum="sre", reference_sets=[self.reference_set]
        )

    def test_panel_generator(self):
        pg = AlleleGenerator(reference_filepath=f"{DATA_DIR}/BX571856.1.fasta", kmer=31)
        assert pg.ref is not None

    def test_simple_snp_variant(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v)
        assert panel.refs[0][:31] != panel.alts[0][:31]
        assert panel.refs[0][-32:] != panel.alts[0][-32:]
        assert panel.refs[0][-31:] != panel.alts[0][-31:]

        assert_no_overlapping_kmers(panel)

        assert panel.refs == [
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
        ]
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTG"
        ]
        assert self.pg._calculate_length_delta_from_indels(v, []) == 0
        assert v.is_indel is False

    def test_simple_variant2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v)
        assert_no_overlapping_kmers(panel)

        assert panel.refs == [
            "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"
        ]
        assert panel.alts == [
            "GATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTGA"
        ]

    def test_simple_variant_invalid(self):
        with pytest.raises(ValueError) as cm:
            v = Variant.create(
                variant_sets=self.variant_sets,
                reference=self.reference,
                reference_bases="T",
                start=31,
                alternate_bases=["T"],
            )
            panel = self.pg.create(v)

    def test_simple_variant_start(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=1,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v)
        #        assert_no_overlapping_kmers(panel) ## Will have overlapping kmers only if the SNP is in the i
        assert panel.refs == [
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
        ]
        assert panel.alts == [
            "TGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG"
        ]

    def test_simple_variant_end(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=2902618,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v)
        assert_no_overlapping_kmers(panel)

        assert panel.refs == [
            "TTTATACTACTGCTCAATTTTTTTACTTTTATNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        ]
        assert panel.alts == [
            "TTTATACTACTGCTCAATTTTTTTACTTTTTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        ]

        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="T",
            start=2902616,
            alternate_bases=["C"],
        )
        panel = self.pg.create(v)
        assert panel.refs == [
            "ATTTTATACTACTGCTCAATTTTTTTACTTTTATNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        ]
        assert panel.alts == [
            "ATTTTATACTACTGCTCAATTTTTTTACTTCTATNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        ]

    def test_simple_variant_with_nearby_snp(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["T"],
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

        assert set(panel.refs) == set(
            [
                "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )
        assert set(panel.alts) == set(
            [
                "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )

    def test_simple_variant_with_multiple_nearby_snps(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["T"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["T"],
        )
        v3 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=30,
            alternate_bases=["G"],
        )

        panel = self.pg.create(v, context=[v2, v3])
        assert_no_overlapping_kmers(panel)

        assert panel.refs == [
            "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTAAAGATAGAAATACACGATGCGAGGAATCAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTAAAGATAGAAATACACGATGCGAGGATTCAAATTTCATAACATCACCATGAGTTTGAT",
        ]
        assert panel.alts == [
            "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTGAT",
        ]

    def test_simple_variant_with_multiple_nearby_snps2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["T"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["T"],
        )
        v3 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=30,
            alternate_bases=["G"],
        )
        v4 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=30,
            alternate_bases=["T"],
        )
        v5 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=30,
            alternate_bases=["A"],
        )

        assert sorted(self.pg._split_context([v, v3, v4])) == sorted([[v, v4], [v, v3]])
        assert (self.pg._split_context([v3, v4])) == [[v4], [v3]]
        assert (self.pg._split_context([v, v3, v4, v5])) == [[v, v4, v5], [v, v3, v5]]
        panel = self.pg.create(v, context=[v2, v3, v4, v5])
        assert_no_overlapping_kmers(panel)
        assert sorted(panel.refs) == sorted(
            [
                "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGATTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTATTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGAAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGAATTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )
        assert sorted(panel.alts) == sorted(
            [
                "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTTTTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGATATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGATTTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )

    def test_simple_variant_with_multiple_nearby_snps(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=31,
            alternate_bases=["T"],
        )
        v2 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["T"],
        )
        v5 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="A",
            start=32,
            alternate_bases=["G"],
        )
        v3 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=30,
            alternate_bases=["G"],
        )
        v4 = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=30,
            alternate_bases=["T"],
        )
        panel = self.pg.create(v, context=[v2, v3, v4, v5])
        assert_no_overlapping_kmers(panel)
        assert sorted(panel.refs) == sorted(
            [
                "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGATTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTAATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTATTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCAGTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGAGTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTAGTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )
        assert sorted(panel.alts) == sorted(
            [
                "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTTATCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTTTTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGCTGTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGGTGTCAAATTTCATAACATCACCATGAGTTTG",
                "CGATTAAAGATAGAAATACACGATGCGAGTTGTCAAATTTCATAACATCACCATGAGTTTG",
            ]
        )
