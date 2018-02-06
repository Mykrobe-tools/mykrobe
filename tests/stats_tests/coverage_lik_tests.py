from mykatlas.stats import percent_coverage_from_expected_coverage
from mykatlas.stats import log_factorial
from mykatlas.stats import log_lik_depth

from math import log
from math import exp

from nose.tools import assert_raises


def test_percentage_coverage():
    assert percent_coverage_from_expected_coverage(
        100) > percent_coverage_from_expected_coverage(10)
    assert percent_coverage_from_expected_coverage(100) == 1
    assert percent_coverage_from_expected_coverage(1) < 1


def test_log_factorial():
    assert log_factorial(4) - log(24) < 0.0001
    assert log_factorial(4) - 3.17 < 0.01
    assert log_factorial(4) - (log(1) + log(2) + log(3) + log(4)) < 0.0001


def test_log_lik_depth():
    assert exp(
        log_lik_depth(
            expected_depth=10,
            depth=10)) > exp(
        log_lik_depth(
            expected_depth=10,
            depth=1))
    assert exp(
        log_lik_depth(
            expected_depth=10,
            depth=10)) > exp(
        log_lik_depth(
            expected_depth=10,
            depth=8))
    assert exp(
        log_lik_depth(
            expected_depth=10,
            depth=10)) == exp(
        log_lik_depth(
            expected_depth=10,
            depth=9))
    assert exp(
        log_lik_depth(
            expected_depth=10,
            depth=10)) > exp(
        log_lik_depth(
            expected_depth=10,
            depth=11))
    assert log_lik_depth(
        expected_depth=100,
        depth=50) < log_lik_depth(
        expected_depth=10,
        depth=9)
    with assert_raises(ValueError) as cm:
        log_lik_depth(expected_depth=0, depth=0)
    with assert_raises(ValueError) as cm:
        log_lik_depth(expected_depth=-1, depth=9)
    with assert_raises(ValueError) as cm:
        log_lik_depth(expected_depth=12, depth=-1)
    with assert_raises(ValueError) as cm:
        log_lik_depth(expected_depth=0, depth=1)
    assert log_lik_depth(expected_depth=1, depth=0) == -1


# TODO. Expect a higher % coverage if k is lower
