import pytest
from mykrobe.typing import ProbeCoverage
from mykrobe.typing import VariantProbeCoverage
from mykrobe.typing import VariantTyper


@pytest.fixture()
def vt():
    return VariantTyper(expected_depths=[100])


@pytest.fixture()
def vt_no_contaim():
    return VariantTyper(expected_depths=[100], contamination_depths=[])


def test_wt_vars(vt):
    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=100, median_depth=100, k_count=100, klen=31
    )
    alternate_coverages = [
        ProbeCoverage(
            min_depth=100, percent_coverage=3, median_depth=100, k_count=3, klen=31
        )
    ]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )
    call = vt.type([v1])
    assert call["genotype"] == [0, 0]
    assert call["info"].get("expected_depths") == [100]


def test_alt_vars(vt):
    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=3, median_depth=100, k_count=3, klen=31
    )
    alternate_coverages = [
        ProbeCoverage(
            min_depth=100, percent_coverage=100, median_depth=100, k_count=100, klen=31
        )
    ]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )
    call = vt.type([v1])
    assert call["genotype"] == [1, 1]


def test_mixed_vars(vt):
    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=100, median_depth=50, k_count=50, klen=31
    )
    alternate_coverages = [
        ProbeCoverage(
            min_depth=100, percent_coverage=100, median_depth=50, k_count=50, klen=31
        )
    ]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )
    call = vt.type(v1)
    assert call["genotype"] == [0, 1]


def test_mixed_vars2(vt):
    reference_coverage = ProbeCoverage(
        min_depth=11, percent_coverage=100, median_depth=42, k_count=42, klen=31
    )
    alternate_coverages = [
        ProbeCoverage(
            min_depth=94, percent_coverage=100, median_depth=102, k_count=94, klen=31
        )
    ]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )
    call = vt.type(v1)
    assert call["genotype"] == [0, 1]


def test_typer_with_contm_simple_case(vt_no_contaim):
    # To do add contamination type
    # self.vt_contaim = VariantTyper(
    #     expected_depths=[80],
    #     contamination_depths=[20])

    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=100, median_depth=80, k_count=80, klen=31
    )
    alternate_coverages = [
        ProbeCoverage(
            min_depth=100, percent_coverage=100, median_depth=20, k_count=40, klen=31
        )
    ]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )

    call = vt_no_contaim.type(v1)
    assert call["genotype"] == [0, 1]

    # call = self.vt_contaim.type(v1)
    # assert call['genotype'] == [0, 0]


def test_typer_with_mult_alt_coverages_simple_case():
    # to do, test should pass on kc model also
    vt_no_contaim = VariantTyper(
        expected_depths=[100], contamination_depths=[], model="median_depth"
    )

    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=70, median_depth=80, k_count=80, klen=31
    )
    alt1 = ProbeCoverage(
        min_depth=100, percent_coverage=70, median_depth=20, k_count=20, klen=31
    )
    alt2 = ProbeCoverage(
        min_depth=100, percent_coverage=100, median_depth=80, k_count=80, klen=31
    )
    alternate_coverages = [alt1, alt2]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )
    assert v1._choose_best_alternate_coverage() == alt2

    call = vt_no_contaim.type(v1)
    assert call["genotype"] == [1, 1]


def test_typer_with_mult_probe_coverages_simple_case(vt_no_contaim):
    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=80, median_depth=80, k_count=80, klen=31
    )
    alt1 = ProbeCoverage(
        min_depth=100, percent_coverage=50, median_depth=20, k_count=20, klen=31
    )
    alt2 = ProbeCoverage(
        min_depth=100, percent_coverage=40, median_depth=80, k_count=30, klen=31
    )
    alternate_coverages = [alt1, alt2]

    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )

    reference_coverage = ProbeCoverage(
        min_depth=100, percent_coverage=80, median_depth=80, k_count=20, klen=31
    )
    alt1 = ProbeCoverage(
        min_depth=100, percent_coverage=50, median_depth=20, k_count=20, klen=31
    )
    alt2 = ProbeCoverage(
        min_depth=100, percent_coverage=100, median_depth=80, k_count=100, klen=31
    )

    alternate_coverages = [alt1, alt2]

    v2 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )

    call = vt_no_contaim.type([v1, v2])
    assert call["genotype"] == [1, 1]


def test_type_with_low_minimum_1(vt_no_contaim):
    reference_coverage = ProbeCoverage(
        min_depth=131, percent_coverage=95.2381, median_depth=155, k_count=131, klen=31
    )
    alt1 = ProbeCoverage(
        min_depth=1, percent_coverage=100, median_depth=1, k_count=1, klen=31
    )
    alternate_coverages = [alt1]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )

    call = vt_no_contaim.type(v1)
    assert call["genotype"] == [0, 0]


def test_type_with_low_minimum_2():
    vt2_no_contaim = VariantTyper(expected_depths=[1], contamination_depths=[])
    reference_coverage = ProbeCoverage(
        min_depth=2, percent_coverage=59.52, median_depth=2, k_count=60, klen=31
    )
    alt1 = ProbeCoverage(
        min_depth=1, percent_coverage=83.33, median_depth=1, k_count=83, klen=31
    )
    alternate_coverages = [alt1]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )

    call = vt2_no_contaim.type(v1)
    assert call["genotype"] == [1, 1]
    assert call["info"]["conf"] < 150


def test_type_with_low_minimum_3():
    vt = VariantTyper(
        expected_depths=[6], contamination_depths=[], confidence_threshold=3
    )
    reference_coverage = ProbeCoverage(
        min_depth=1, percent_coverage=100, median_depth=2, k_count=2, klen=31
    )
    alt1 = ProbeCoverage(
        min_depth=1, percent_coverage=100, median_depth=1, k_count=1, klen=31
    )
    alternate_coverages = [alt1]
    v1 = VariantProbeCoverage(
        var_name="A123T",
        reference_coverages=[reference_coverage],
        alternate_coverages=alternate_coverages,
    )

    call = vt.type(v1)
    assert call["genotype"] == [0, 1]
    assert call["info"]["conf"] < 100
