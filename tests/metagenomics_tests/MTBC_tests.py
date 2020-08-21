from unittest import TestCase
import os
import sys
sys.path.append(".")
from mykrobe.metagenomics import AMRSpeciesPredictor
DATA_DIR = os.path.join("tests", "ref_data")

class MTBCSpeciesTests(TestCase):

    def setUp(self):
        self.hierarchy_json_file = f"{DATA_DIR}/mtbc_hierarchy.json"

    def teardown(self):
        pass

    def test_mixed_MTBC_NTM(self):
        species_predictor = AMRSpeciesPredictor(
            phylo_group_covgs={},
            sub_complex_covgs={},
            species_covgs={},
            lineage_covgs={},
            hierarchy_json_file=self.hierarchy_json_file)
        species_predictor.out_json["phylogenetics"] = {}
        species_predictor.out_json["phylogenetics"]["phylo_group"] = {
            "Non_tuberculosis_mycobacterium_complex": {
                "percent_coverage": 58.71542975006994,
                "median_depth": 36
            },
            "Mycobacterium_tuberculosis_complex": {
                "percent_coverage": 62.81850563578579,
                "median_depth": 2
            }
        }
        assert species_predictor.is_mtbc_present()
        assert species_predictor.is_ntm_present()

        assert len(species_predictor._get_present_phylo_groups(
            species_predictor.out_json["phylogenetics"]["phylo_group"])) == 2

    def test_get_best_coverage_dict(self):
        species_predictor = AMRSpeciesPredictor(
            phylo_group_covgs={},
            sub_complex_covgs={},
            species_covgs={},
            lineage_covgs={},
            hierarchy_json_file=self.hierarchy_json_file)

        best_species = species_predictor._get_best_coverage_dict({
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
            }}).keys()
        assert list(best_species) == ["Mycobacterium_chimaera"]

    def test_mixed_chimera(self):
        species_predictor = AMRSpeciesPredictor(
            phylo_group_covgs={},
            sub_complex_covgs={},
            species_covgs={},
            lineage_covgs={},
            hierarchy_json_file=self.hierarchy_json_file)
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
