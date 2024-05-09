from base import assert_no_overlapping_kmers
from mykrobe.variants.schema.models import Variant


def test_simple_deletion1(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AA",
        start=31,
        alternate_bases=["A"],
    )
    assert v.is_indel
    assert v.is_deletion
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert pg._calculate_length_delta_from_indels(v, []) == 1
    assert panel.alts == ["GATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTG"]


def test_simple_deletion2(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AT",
        start=32,
        alternate_bases=["A"],
    )
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert panel.alts == [
        "ATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGAT"
    ]


def test_simple_deletion3(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AT",
        start=2902618,
        alternate_bases=["T"],
    )
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "TTTATACTACTGCTCAATTTTTTTACTTTTATNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" in panel.refs
    assert panel.alts == [
        "TTTATACTACTGCTCAATTTTTTTACTTTTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    ]


def test_simple_deletion4(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="ATC",
        start=32,
        alternate_bases=["A"],
    )
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert panel.alts == ["ATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGAT"]


def test_simple_insertion1(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="C",
        start=1,
        alternate_bases=["TTTC"],
    )
    panel = pg.create(v)
    #        assert_no_overlapping_kmers(panel)### Skip this test for vars in first k bases of ref
    assert v.is_indel
    assert v.is_insertion
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert panel.alts == ["TTTCGATTAAAGATAGAAATACACGATGCGAGC"]


def test_simple_insertion2(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="C",
        start=1,
        alternate_bases=["CTTT"],
    )
    panel = pg.create(v)
    #        assert_no_overlapping_kmers(panel)### Skip this test for vars in first k bases of ref
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert panel.alts == ["CTTTGATTAAAGATAGAAATACACGATGCGAGCA"]


def test_simple_insertion3(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=31,
        alternate_bases=["ATTT"],
    )
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTG" in panel.refs
    assert panel.alts == [
        "GATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTG"
    ]


def test_simple_insertion4(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=32,
        alternate_bases=["AGGGG"],
    )
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA" in panel.refs
    assert panel.alts == [
        "ATTAAAGATAGAAATACACGATGCGAGCAAGGGGTCAAATTTCATAACATCACCATGAGTTTGA"
    ]


def test_simple_insertion5(variant_sets_and_reference, pg):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=2902618,
        alternate_bases=["ATGC"],
    )
    panel = pg.create(v)
    assert_no_overlapping_kmers(panel)
    assert "TTTATACTACTGCTCAATTTTTTTACTTTTATNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" in panel.refs
    assert panel.alts == [
        "TATACTACTGCTCAATTTTTTTACTTTTATGCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    ]


def test_double_insertion(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=4021408,
        alternate_bases=["ACGCTGGCGGGCG"],
    )
    v1 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AGA",
        start=4021406,
        alternate_bases=["CGG"],
    )
    context = [v1]
    assert pg2._remove_overlapping_contexts(v, [v1]) == []
    panel = pg2.create(v, context=context)
    assert_no_overlapping_kmers(panel)
    assert "ATCTAGCCGCAAGGGCGCGAGCAGACGCAGAATCGCATGATTTGAGCTCAAATCATGCGAT" in panel.refs
    assert panel.alts == [
        "TCTAGCCGCAAGGGCGCGAGCAGACGCAGACGCTGGCGGGCGATCGCATGATTTGAGCTCAAATCATGCGAT"
    ]


def test_double_indel_fail(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="CCA",
        start=2288851,
        alternate_bases=["A"],
    )
    v1 = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="A",
        start=2288850,
        alternate_bases=["ACC"],
    )
    context = [v1]
    panel = pg2.create(v, context=context)
    assert "GGCGCACACAATGATCGGTGGCAATACCGACCACATCGACCTCATCGACGCCGCGTTGCCG" in panel.refs
    assert (
        "GGCGCACACAATGATCGGTGGCAATACCGACCACATCGACCTCATCGACGCCGCGTTGCCG"
        not in panel.alts
    )


def test_large_insertion(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="CCGCCGGCCCCGCCGTTT",
        start=1636155,
        alternate_bases=[
            "CTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCG"
        ],
    )
    panel = pg2.create(v, context=[])
    assert_no_overlapping_kmers(panel)
    assert "AGACCTAGCAGGGTGCCGGCGCCGCCCTTGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCAT" in panel.refs
    assert panel.alts == [
        "GACCTAGCAGGGTGCCGGCGCCGCCCTTGCTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCGCCATCGCCGATGATGTTTTCC"
    ]
