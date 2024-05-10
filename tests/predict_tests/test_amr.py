import os

from mykrobe.predict import TBPredictor
from mykrobe.variants.schema.models import Variant


def test_wr_vars():
    variant_snp = Variant.create(
        start=0, end=1, reference_bases="A", alternate_bases=["T"]
    )

    predictor = TBPredictor(
        variant_calls={},
        called_genes={},
        variant_to_resistance_json_fp=os.path.join(
            "tests", "ref_data", "tb_variant_to_resistance_drug.json"
        ),
    )

    call = {
        "variant": None,
        "genotype": [0, 1],
        "genotype_likelihoods": [0.1, 0.9, 0.12],
        "info": {
            "contamination_depths": [],
            "coverage": {
                "alternate": {
                    "percent_coverage": 100.0,
                    "median_depth": 15,
                    "min_depth": 2,
                },
                "reference": {
                    "percent_coverage": 100.0,
                    "median_depth": 139,
                    "min_depth": 128,
                },
            },
            "expected_depths": [152],
        },
    }

    assert predictor._coverage_greater_than_threshold(call, [""]) is False
