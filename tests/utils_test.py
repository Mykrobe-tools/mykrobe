import pytest

from mykrobe import utils


def test_x_mutation_fixed_var_name():
    assert "katG_S315G-GCT2155167GGT" == utils._x_mutation_fixed_var_name(
        "katG_S315X-GCT2155167GGT"
    )
    assert utils._x_mutation_fixed_var_name("katG_S315.-GCT2155167GGT") is None
    assert utils._x_mutation_fixed_var_name("katG_S315X-GCT2155167NNN") is None


def test_fix_amino_acid_X_variants_keys():
    test_dict = {
        "foo": "bar",
        "katG_S315X-GCT2155167GGT": "baz",
        "katG_S315C-GCT2155167CTT": "baz",
    }

    utils.fix_amino_acid_X_variants_keys(test_dict)
    assert test_dict == {
        "foo": "bar",
        "katG_S315G-GCT2155167GGT": "baz",
        "katG_S315C-GCT2155167CTT": "baz",
    }
