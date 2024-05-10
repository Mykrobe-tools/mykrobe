import pytest

from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantSet
from mykrobe.variants.schema.models import VariantCall
from mykrobe.variants.schema.models import VariantCallSet


@pytest.fixture()
def call_set(reference_set, reference):
    variant_set = VariantSet.create_and_save(
        name="this_vcf_file2", reference_set=reference_set
    )
    variant_sets = [variant_set]
    return VariantCallSet.create(
        sample_id="C00123", name="C00123", variant_sets=variant_sets
    )


@pytest.fixture()
def variant_snp(variant_sets, reference):
    return Variant.create(
        variant_sets=variant_sets,
        start=0,
        end=1,
        reference_bases="A",
        alternate_bases=["T"],
        reference=reference,
    )


def test_create_call_set(variant_sets):
    call_set = VariantCallSet.create_and_save(
        name="call_set", sample_id="C00123", variant_sets=variant_sets
    )
    cs = VariantCallSet.objects.get(name="call_set")
    assert call_set == cs
    assert cs.name == "call_set"
    assert cs.variant_sets[0].reference_set.name == "ref_set"


def test_create_SNP_het_call(variant_snp, call_set):
    c1 = VariantCall.create(
        variant=variant_snp,
        call_set=call_set,
        genotype=[0, 1],
        genotype_likelihoods=[0.1, 0.9, 0.12],
    )
    assert c1.call_set_name == "C00123"
    assert c1.genotype == [0, 1]
    variant_snp.save()
    c1.save()
    assert c1 in variant_snp.calls

    c2 = VariantCall.create(
        variant=variant_snp,
        call_set=call_set,
        genotype="1/1",
        genotype_likelihoods=[0.01, 0.1, 0.9],
    )
    assert c2.call_set_name == "C00123"
    assert c2.genotype == [1, 1]

    c2.save()
    assert c2 in variant_snp.calls


def test_create_complex_call(reference, variant_snp, variant_sets, call_set):
    variant_snp_mult_alts = Variant.create(
        variant_sets=variant_sets,
        start=0,
        end=1,
        reference_bases="T",
        alternate_bases=["A", "C"],
        reference=reference,
    )
    c1 = VariantCall.create(
        variant=variant_snp_mult_alts,
        call_set=call_set,
        genotype="2/1",
        genotype_likelihoods=[0.01, 0.1, 0.9, 0.1, 0.2, 0.6],
    )
    variant_snp.save()
    variant_snp_mult_alts.save()
    c1.save()
    assert c1.call_set_name == "C00123"
    assert c1.genotype == [2, 1]
    assert c1 not in variant_snp.calls
    assert c1 in variant_snp_mult_alts.calls
