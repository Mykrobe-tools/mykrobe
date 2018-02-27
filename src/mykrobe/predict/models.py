# Defines Mykrobe predictor output
from mongoengine import Document
from mongoengine import DictField

from mykrobe.utils import unique


class MykrobePredictorSusceptibilityResult(Document):

    susceptibility = DictField()

    @classmethod
    def create(cls, susceptibility):
        return cls(susceptibility=susceptibility)

    def to_dict(self):
        return self.to_mongo().to_dict()

    def __eq__(self, other):
        return self.susceptibility == other.susceptibility

    def diff(self, other):
        diff = {}
        # Compares the antibiogram of two predictor results
        drugs = unique(self.drugs + other.drugs)
        for drug in drugs:
            predict1, predict2 = self.susceptibility.get(drug, {"predict": "NA"}).get(
                "predict"), other.susceptibility.get(drug, {"predict": "NA"}).get("predict")
            if predict1 != predict2:
                diff[drug] = {}
                diff[drug]["predict"] = (predict1, predict2)
        return diff

    @property
    def drugs(self):
        return list(self.susceptibility.keys())
