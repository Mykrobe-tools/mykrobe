from unittest import TestCase
from mykrobe.metagenomics import AMRSpeciesPredictor


class MTBCSpeciesTests(TestCase):

    def setUp(self):
        pass

    def teardown(self):
        pass

    def test_mixed_chimera(self):
        species_predictor = AMRSpeciesPredictor(
            phylo_group_covgs={},
            sub_complex_covgs={},
            species_covgs={},
            lineage_covgs={}
        )
        species_predictor.out_json["phylogenetics"] = {
            "sub_complex": {
                "Mycobacterium_avium_complex": {
                    "percent_coverage": 98.346,
                    "median_depth": 54.0
                }
            },
            "phylo_group": {
                "Non_tuberculosis_mycobacterium_complex": {
                    "percent_coverage": 82.846,
                    "median_depth": 49
                }
            },
            "species": {
                "Mycobacterium_chimaera": {
                    "percent_coverage": 99.162,
                    "median_depth": 39
                },
                "Mycobacterium_intracellulare": {
                    "percent_coverage": 98.662,
                    "median_depth": 45
                },
                "Mycobacterium_bovis": {
                    "percent_coverage": 9.894,
                    "median_depth": 12.0
                }
            }
        }

        out_dict = species_predictor.choose_best(
            species_predictor.out_json["phylogenetics"])

        assert "Mycobacterium_chimaera" in out_dict["species"]
        assert "Mycobacterium_intracellulare" in out_dict["species"]
        assert "Mycobacterium_bovis" not in out_dict["species"]
