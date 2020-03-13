import sys
import pytest

sys.path.append(".")
from unittest import TestCase

from mykrobe.variants.schema.models import VariantCall
from mykrobe.variants.schema.models import Variant
from mykrobe.cmds.amr import Panel, Species, TbPanel, StaphPanel, GnPanel, PANELS

from mykrobe.predict import TBPredictor


class AMRPredictTest(TestCase):
    def setUp(self):
        self.variant_snp = Variant.create(
            start=0, end=1, reference_bases="A", alternate_bases=["T"]
        )

        self.predictor = TBPredictor(variant_calls={}, called_genes={})

    def teardown(self):
        pass

    def test_wt_vars(self):
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

        assert self.predictor._coverage_greater_than_threshold(call, [""]) == False


class TestPanel:
    def test_fromSpeciesAndName_invalidSpecies_raisesError(self):
        species = "foo"
        name = "bar"
        with pytest.raises(NameError):
            Panel.from_species_and_name(species, name)

    def test_fromSpeciesAndName_tbWithInvalidPanel_raisesError(self):
        species = Species.TB
        name = "bar"
        with pytest.raises(ValueError):
            Panel.from_species_and_name(species, name)

    def test_fromSpeciesAndName_staphWithInvalidPanel_raisesError(self):
        species = Species.STAPH
        name = "bar"
        with pytest.raises(ValueError):
            Panel.from_species_and_name(species, name)

    def test_fromSpeciesAndName_gnWithInvalidPanel_raisesError(self):
        species = Species.GN
        name = "bar"
        with pytest.raises(ValueError):
            Panel.from_species_and_name(species, name)

    def test_fromSpeciesAndName_staphWithCustomPanel_returnsCustom(self):
        species = Species.STAPH
        name = "custom"
        panel = Panel.from_species_and_name(species, name)

        actual = panel.name
        expected = StaphPanel.CUSTOM

        assert actual == expected

    def test_fromSpeciesAndName_TbWithValidPanel_returnsPanel(self):
        species = Species.TB
        name = "201901"
        panel = Panel.from_species_and_name(species, name)

        actual = panel.name
        expected = TbPanel.NEJM_WALKER

        assert actual == expected

    def test_fromSpeciesAndName_GnWithValidPanel_returnsPanel(self):
        species = Species.GN
        name = "default"
        panel = Panel.from_species_and_name(species, name)

        actual = panel.name
        expected = GnPanel.DEFAULT

        assert actual == expected

    def test_add_multiple_paths(self):
        species = Species.TB
        name = "201901"
        panel = Panel.from_species_and_name(species, name)
        additions = ["foo", "bar"]
        panel.add_path(*additions)

        actual = panel.paths
        expected = PANELS[species][panel.name]
        expected.extend(additions)

        assert actual == expected

    def test_add_single_path(self):
        species = Species.TB
        name = "201901"
        panel = Panel.from_species_and_name(species, name)
        additions = "foo"
        panel.add_path(*additions)

        actual = panel.paths
        expected = PANELS[species][panel.name]
        expected.append(additions)

        assert actual == expected
