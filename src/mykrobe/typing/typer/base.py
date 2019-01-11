MIN_LLK = -99999999
DEFAULT_ERROR_RATE = 0.05
DEFAULT_MINOR_FREQ = 0.2
import logging
logger = logging.getLogger(__name__)


class Typer(object):

    def __init__(
            self,
            expected_depths,
            contamination_depths=[],
            error_rate=0.05,
            ignore_filtered=False,
            filters=[],
            confidence_threshold=1
    ):
        self.expected_depths = expected_depths
        self.contamination_depths = contamination_depths
        self.error_rate = error_rate
        self.ignore_filtered = ignore_filtered
        self.filters = filters
        self.confidence_threshold = confidence_threshold

    def type(self, l):
        raise NotImplemented("Implemented in sub class")

    def likelihoods_to_genotype(self, likelihoods, min_like=MIN_LLK):
        ml = max(likelihoods)
        i = likelihoods.index(ml)
        if i == 0:
            if ml <= min_like:
                gt = "-/-"
            else:
                gt = "0/0"
        elif i == 1:
            gt = "0/1"
        elif i == 2:
            gt = "1/1"
        return gt

    def has_contamination(self):
        return self.contamination_depths or len(self.expected_depths) > 1
