from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantSet
from mykrobe.utils import split_var_name


def test_create_new_variant(reference_set):
    variant_set = VariantSet.create_and_save(
        name="this_vcf_file", reference_set=reference_set
    )
    vs = VariantSet.objects.get(name="this_vcf_file")
    assert variant_set == vs
    assert vs.reference_set.name == "ref_set"


def test_create_SNP(reference, variant_sets):
    v1 = Variant.create(
        variant_sets=variant_sets,
        start=0,
        end=1,
        reference_bases="A",
        alternate_bases=["T"],
        reference=reference,
    )
    assert v1.start == 0
    assert v1.end == 1
    assert v1.alternate_bases == ["T"]
    assert v1.length == 0


def test_create_insertion(reference, variant_sets):
    v1 = Variant.create(
        variant_sets=variant_sets,
        start=0,
        end=1,
        reference_bases="T",
        alternate_bases=["TA"],
        reference=reference,
    )
    assert v1.start == 0
    assert v1.end == 1
    assert v1.alternate_bases == ["TA"]
    assert v1.is_insertion
    assert v1.is_deletion is False
    assert v1.is_indel
    assert v1.length == 1


def test_create_deletion(reference, variant_sets):
    v1 = Variant.create(
        variant_sets=variant_sets,
        start=0,
        end=1,
        reference_bases="AA",
        alternate_bases=["A"],
        reference=reference,
    )
    assert v1.start == 0
    assert v1.end == 1
    assert v1.alternate_bases == ["A"]
    assert v1.reference_bases == "AA"
    assert v1.is_insertion is False
    assert v1.is_deletion
    assert v1.is_indel
    assert v1.length == 1


def test_split_name():
    name = "A12T"
    r, pos, a = split_var_name(name)
    assert r == "A"
    assert pos == 12
    assert a == "T"


def test_split_name_del():
    name = "AA12T"
    r, pos, a = split_var_name(name)
    assert r == "AA"
    assert pos == 12
    assert a == "T"


def test_split_name_ins():
    name = "A12TT"
    r, pos, a = split_var_name(name)
    assert r == "A"
    assert pos == 12
    assert a == "TT"


def test_split_name2():
    name = "A12T/A"
    r, pos, a = split_var_name(name)
    assert r == "A"
    assert pos == 12
    assert a == "T/A"


def test_split_name3():
    name = "C-54T"
    r, pos, a = split_var_name(name)
    assert r == "C"
    assert pos == -54
    assert a == "T"
