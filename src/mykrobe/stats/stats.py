from math import exp
from math import factorial
from math import log
import logging
logger = logging.getLogger(__name__)


def percent_coverage_from_expected_coverage(coverage):
    # With low coverage we expect a lower percent of the sequence to be
    # coverage.
    return 1 - exp(-coverage)


def log_lik_probability_of_N_gaps(depth, percent_coverage, L):
    percent_coverage = float(percent_coverage)/100
    n_gaps = int(round(L-(L*percent_coverage)))
    expected_n_gaps = exp(-depth) * L
    try:
        log_like= log_poisson_prob(expected_n_gaps, n_gaps)
    except ValueError:
        ## For very high coverage samples, expected_n_gaps goes to 0.0 breaking math.log
        if expected_n_gaps<=0:
            expected_n_gaps=1e-308        
        log_like= log_poisson_prob(expected_n_gaps, n_gaps)
    return log_like


def log_poisson_prob(lam, k):
    return -lam + k * log(lam) - log_factorial(k)


def log_factorial(n):
    assert n >= 0
    out = 0
    for i in range(int(n)):
        out += log(i + 1)
    return out


def log_lik_depth(depth, expected_depth):
    if expected_depth <= 0:
        raise ValueError("Expected depth must be greater than 0")
    if depth < 0:
        raise ValueError("Depth must not be negative")
    return log_poisson_prob(lam=expected_depth, k=depth)


def log_lik_R_S_coverage(observed_alternate_depth,
                         observed_reference_depth,
                         expected_alternate_depth,
                         expected_reference_depth):
    lne = log_poisson_prob(
        lam=expected_alternate_depth,
        k=observed_alternate_depth)
    le = log_poisson_prob(
        lam=expected_reference_depth,
        k=observed_reference_depth)
    return lne + le


def log_lik_R_S_kmer_count(observed_reference_kmer_count,
                           observed_alternate_kmer_count,
                           expected_reference_kmer_count,
                           expected_alternate_kmer_count):
    # logger.debug("%f, %f, %f" % (expected_reference_depth,
    #                              expected_reference_kmer_count, observed_reference_kmer_count))
    # logger.debug("%f, %f" % (expected_reference_kmer_count,expected_alternate_kmer_count))
    # logger.debug("%f, %f" % (observed_reference_kmer_count,observed_alternate_kmer_count))
    lne = log_poisson_prob(
        lam=expected_reference_kmer_count, k=observed_reference_kmer_count)
    le = log_poisson_prob(
        lam=expected_alternate_kmer_count, k=observed_alternate_kmer_count)
    # logger.debug("%i, %i, %i, %f" % (expected_reference_depth,
    #                                  expected_reference_kmer_count, observed_reference_kmer_count, lne))
    # logger.debug("%i, %i, %i, %f" % (expected_alternate_depth,
    # expected_alternate_kmer_count, observed_alternate_kmer_count, le))
    return lne + le
