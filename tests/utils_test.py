from io import StringIO

import pytest
from Bio.Seq import Seq

from mykrobe import utils
from mykrobe.utils import get_first_chrom_name, is_file_compressed, load_json


# When we have an amino acid to any other amino mutation (eg I42X), we want to
# replace the X with the actual amino acid. This means getting the codon from
# the variant name and translating it. We have the variants in the form eg:
#   katG_S315X-GCT2155167GGT
# Problem is, if the gene is on the reverse strand then we have to reverse
# complement the codon. Now, we don't actually know the strand that the gene is
# on, but we do know the first codon and the amino acid - in that example they
# are GCT and S. This test is not testing the mykrobe code, but is testing
# that we can use the first codon to check if a gene is on the reverse strand.
# For this to work, we can't have a codon and its reverse complement translating
# to the same amino acid.
# Eg for GCT -> S, the reverse complement is AGC, which translates to A. This
# test runs through all codons and checks that we never get the same amino
# acid from a codon and the codon's reverse complement.


def test_x_mutation_no_amino_acid_revcomp_clash():
    nucs = ["A", "C", "G", "T"]
    codons = [f"{x}{y}{z}" for x in nucs for y in nucs for z in nucs]

    for codon in codons:
        amino_acid = str(Seq(codon).translate())
        amino_acid_rev = str(Seq(codon).reverse_complement().translate())
        assert amino_acid != amino_acid_rev


def test_x_mutation_fixed_var_name():
    # This is on the reverse strand, so the GGT should get reverse complemented,
    # resulting in amino acid T (=ACC), not G (=GGT).
    assert "katG_S315T-GCT2155167GGT" == utils._x_mutation_fixed_var_name(
        "katG_S315X-GCT2155167GGT"
    )
    assert utils._x_mutation_fixed_var_name("katG_S315.-GCT2155167GGT") is None
    assert utils._x_mutation_fixed_var_name("katG_S315X-GCT2155167NNN") is None

    # This is on forward strand, so should get TAT -> Y
    assert "embB_D328Y-GAT4247495TAT" == utils._x_mutation_fixed_var_name(
        "embB_D328X-GAT4247495TAT"
    )


def test_fix_amino_acid_X_variants_keys():
    test_dict = {
        "foo": "bar",
        "katG_S315X-GCT2155167GGT": "baz",
        "katG_S315C-GCT2155167TTA": "baz",
        "katG_S315X-GCT2155167CTA": "baz",  # katG is on reverse strand so stop TAG is CTA
        "embB_M306I-ATG4247429ATA": "value",
        "embB_M306X-ATG4247429ATA": "value",  # duplicate of previous, should remove
    }

    utils.fix_amino_acid_X_variants_keys(test_dict)
    assert test_dict == {
        "foo": "bar",
        "katG_S315T-GCT2155167GGT": "baz",
        "katG_S315C-GCT2155167TTA": "baz",
        "katG_S315*-GCT2155167CTA": "baz",
        "embB_M306I-ATG4247429ATA": "value",
    }


def test_get_first_chrom_name_not_fasta_raises_error():
    fp = StringIO("@chrom1\nACGT")
    with pytest.raises(ValueError):
        get_first_chrom_name(fp)


def test_get_first_chrom_name_multiple_chroms_gets_first():
    fp = StringIO(">ch1\nA\n>ch2\nG")

    actual = get_first_chrom_name(fp)
    expected = "ch1"

    assert actual == expected


def test_get_first_chrom_name_pointer_is_past_start_gets_first():
    fp = StringIO(">ch1\nA\n>ch2\nG")
    _ = fp.readline()

    actual = get_first_chrom_name(fp)
    expected = "ch1"

    assert actual == expected


def test_is_compressed_file_file_is_not_compressed():
    filepath = "tests/vcf_tests/test.vcf"
    assert not is_file_compressed(filepath)


def test_is_compressed_file_file_is_compressed():
    filepath = "tests/species_data_tests/species1_data.20200101.tar.gz"
    assert is_file_compressed(filepath)


def test_load_json_uncompressed():
    filepath = "tests/ref_data/tb_variant_to_resistance_drug.json"
    data = load_json(filepath)

    assert data["katG_S315T"] == ["Isoniazid"]


def test_load_json_compressed():
    filepath = "tests/ref_data/tb_variant_to_resistance_drug.json.gz"
    data = load_json(filepath)

    assert data["katG_S315T"] == ["Isoniazid"]
