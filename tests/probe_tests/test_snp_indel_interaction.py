from base import assert_no_overlapping_kmers
from mykrobe.variants.schema.models import Variant


def test_ins_with_SNP_context(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=31,
        alternate_bases=["ATTT"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=32,
        alternate_bases=["T"],
    )
    panel = pg.create(v, context=[v2])
    # assert_no_overlapping_kmers(panel)  ### This test seems to fail sometimes, and pass othertimes...
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "GATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTG",
            "TTAAAGATAGAAATACACGATGCGAGCATTTTTCAAATTTCATAACATCACCATGAGTTTG",
        ]
    )


def test_del_with_SNP_context1(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AA",
        start=31,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=33,
        alternate_bases=["A"],
    )
    panel = pg.create(v, context=[v2])
    assert_no_overlapping_kmers(panel)
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "ATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGA",
            "GATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTG",
        ]
    )


def test_del_with_SNP_context2(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AA",
        start=31,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=32,
        alternate_bases=["T"],
    )
    panel = pg.create(v, context=[v2])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2]) == []
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert sorted(panel.alts) == sorted(
        ["GATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTG"]
    )


def test_del_with_ins_context1(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AAT",
        start=31,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=4,
        alternate_bases=["TTTT"],
    )
    panel = pg.create(v, context=[v2])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2]) == [v2]
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "GATTAAAGATAGAAATACACGATGCGAGCACAAATTTCATAACATCACCATGAGTTTGAT",
            "TTTTAAAGATAGAAATACACGATGCGAGCACAAATTTCATAACATCACCATGAGTTTG",
        ]
    )


def test_del_with_ins_context2(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="C",
        start=1,
        alternate_bases=["CTTT"],
    )
    panel = pg.create(v, context=[v2])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2]) == [v2]
    assert pg._remove_contexts_not_within_k(v, [v2]) == []
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        ["ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT"]
    )


def test_del_with_ins_context3(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=5,
        alternate_bases=["TT"],
    )
    panel = pg.create(v, context=[v2])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2]) == [v2]
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "TTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
        ]
    )


def test_del_with_ins_context4(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=5,
        alternate_bases=["TT"],
    )
    v3 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=5,
        alternate_bases=["TG"],
    )
    panel = pg.create(v, context=[v2, v3])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2, v3]) == [v2, v3]
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "TTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "TTGAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
        ]
    )


def test_del_with_ins_context5(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=5,
        alternate_bases=["TT"],
    )
    v3 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=6,
        alternate_bases=["AG"],
    )
    panel = pg.create(v, context=[v2, v3])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2, v3]) == [v2, v3]
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "TTAGAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGA",
            "TTAGAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "TTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
        ]
    )


def test_del_with_ins_context_where_base_is_deleted1(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=33,
        alternate_bases=["C"],
    )
    panel = pg.create(v, context=[v2])
    assert_no_overlapping_kmers(panel)
    assert pg._remove_overlapping_contexts(v, [v2]) == []
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        ["ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT"]
    )


def test_del_with_ins_context_where_base_is_deleted2(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="TAAA",
        start=5,
        alternate_bases=["T"],
    )
    v3 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=7,
        alternate_bases=["AG"],
    )
    panel = pg.create(v, context=[v2, v3])
    assert_no_overlapping_kmers(panel)
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATC",
            "TTAAGAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
        ]
    )

    panel = pg.create(v, context=[v3, v2])
    assert_no_overlapping_kmers(panel)
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert sorted(panel.alts) == sorted(
        [
            "ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
            "CGATTGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATC",
            "TTAAGAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT",
        ]
    )


def test_snp_with_replace_context(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="G",
        start=2338961,
        alternate_bases=["A"],
    )
    v1 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="GGATG",
        start=2338990,
        alternate_bases=["CGATA"],
    )
    panel = pg2.create(v, context=[v1])
    assert_no_overlapping_kmers(panel)
    assert "CGACTAGCCACCATCGCGCATCAGTGCGAGGTCAAAAGCGACCAAAGCGAGCAAGTCGCGG" in panel.refs

    assert set(panel.alts) == set(
        [
            "CGACTAGCCACCATCGCGCATCAGTGCGAGATCAAAAGCGACCAAAGCGAGCAAGTCGCCG",
            "CGACTAGCCACCATCGCGCATCAGTGCGAGATCAAAAGCGACCAAAGCGAGCAAGTCGCGG",
        ]
    )


def test_indel_snp_indel_context(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="TCGCGTGGC",
        start=4021459,
        alternate_bases=["GCGAGCAGA"],
    )
    v1 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=4021455,
        alternate_bases=["ATCTAGCCGCAAG"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="T",
        start=4021489,
        alternate_bases=["G"],
    )
    panel = pg2.create(v)  # , context = [v1, v2])
    assert_no_overlapping_kmers(panel)
    assert "ATCATGCGATTCTGCGTCTGCTCGCGAGGCTCGCGTGGCCGCCGGCGCTGGCGGGCGATCT" in panel.refs

    panel = pg2.create(v, context=[v1, v2])
    assert_no_overlapping_kmers(panel)
    assert sorted(panel.alts) == sorted(
        [
            "ATCATGCGATTCTGCGTCTGCTCGCGAGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCG",
            "ATCATGCGATTCTGCGTCTGCTCGCGAGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCT",
            "TGCGTCTGCTCGCGATCTAGCCGCAAGGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCG",
            "TGCGTCTGCTCGCGATCTAGCCGCAAGGGCGCGAGCAGACGCCGGCGCTGGCGGGCGATCT",
        ]
    )


def test_complex_context(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATTT",
        start=1503643,
        alternate_bases=["A"],
    )
    v1 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="CCT",
        start=1503615,
        alternate_bases=["C"],
    )
    v2 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=1503655,
        alternate_bases=["ATGCCGCCGCC"],
    )
    panel = pg2.create(v, context=[v1, v2])
    assert_no_overlapping_kmers(panel)
    assert "ATCCTGGAGCCCACCAGCGGAAACACCGGCATTTCGCTGGCGATGGCGGCCCGGTTGAAGG" in panel.refs
    assert set(panel.alts) == set(
        [
            "CCATCGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGGCGGCCCGGTTGAAGGGGT",
            "TCCTGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGGCGGCCCGGTTGAAGGGG",
            "ATCGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGCCGCCGCCTGGCGGCCCGG",
            "TCCTGGAGCCCACCAGCGGAAACACCGGCACGCTGGCGATGCCGCCGCCTGGCGGCCCGG",
        ]
    )
