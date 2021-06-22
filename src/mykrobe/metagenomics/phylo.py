from __future__ import print_function
import json
import os
import operator
from mykrobe.utils import median
from mykrobe.utils import load_json
from mykrobe.utils import flatten
from mykrobe.metagenomics.models import MykrobePredictorPhylogeneticsResult
from mykrobe.stats import percent_coverage_from_expected_coverage


DEFAULT_THRESHOLD = 30

TAXON_THRESHOLDS = {
    "Saureus" : 30,
    "Sepidermidis" : 30,
    "Shaemolyticus" : 30,
    "Sother" : 15,
    "Coagneg" : 30,
    "Staphaureus" : 30,
    "Escherichia_coli" : 15,
    "Klebsiella_pneumoniae" : 15,
    "Ecoli_Shigella": 90,
    "Shigella_sonnei": 90,
    "Salmonella_enterica": 90,
    "Salmonella_Typhi": 90,
}

class Hierarchy(object):

    def __init__(self, _dict):
        self.dict = _dict

    def get_children(self, target_species):
        phylo_group = self.get_phylo_group(target_species)
        return phylo_group.get("children", {})

    def get_phylo_group(self, target_species):
        phylo_group = []
        for k, v in self.dict.items():
            if k == target_species:
                return v
            for k2, v2 in v["children"].items():
                if k2 == target_species:
                    return v2
                for k3, v3 in v2["children"].items():
                    if k3 == target_species:
                        return v3
                    for k4, v4 in v3["children"].items():
                        if k4 == target_species:
                            return v4


class SpeciesPredictor(object):

    def __init__(
            self,
            phylo_group_covgs,
            sub_complex_covgs,
            species_covgs,
            lineage_covgs,
            verbose=False,
            hierarchy_json_file=None):
        self.phylo_group_covgs = phylo_group_covgs
        self.sub_complex_covgs = sub_complex_covgs
        self.species_covgs = species_covgs
        self.lineage_covgs = lineage_covgs
        self.out_json = {}
        self.threshold = {}
        self.verbose = verbose
        try:
            self.hierarchy = Hierarchy(load_json(hierarchy_json_file))
        except TypeError:
            self.hierarchy = {}

    def run(self):
        self._load_taxon_thresholds()
        self._aggregate_all()
        return MykrobePredictorPhylogeneticsResult(phylogenetics=self.out_json["phylogenetics"])

    def _add_unknown_where_empty(self, covgs):
        if not covgs:
            covgs["Unknown"] = {"percent_coverage": -1, "median_depth": -1}

    def _load_taxon_thresholds(self):
        #taxon_coverage_threshold_file = os.path.realpath(
        #    os.path.join(
        #        os.path.dirname(__file__),
        #        "..",
        #        "data/predict/taxon_coverage_threshold.json"))
        #with open(taxon_coverage_threshold_file, "r") as infile:
        #    self.threshold = json.load(infile)
        self.threshold = TAXON_THRESHOLDS

    def calc_expected_depth(self):
        # Get all of the panels with % coverage > 30
        _median = []
        for phylo_group, coverage_dict in self.phylo_group_covgs.items():
            _median.extend(coverage_dict["median"])
        if _median:
            return median(_median)
        else:
            return 0

    def _aggregate_all(self):
        # Calculate expected coverage
        self.expected_depth = self.calc_expected_depth()
        self._aggregate(self.phylo_group_covgs)
        self._aggregate(self.sub_complex_covgs, threshold=50)
        self._aggregate(self.species_covgs)
        self._aggregate(self.lineage_covgs)
        self.out_json["phylogenetics"] = {}
        self.out_json["phylogenetics"]["phylo_group"] = self.phylo_group_covgs
        self.out_json["phylogenetics"]["sub_complex"] = self.sub_complex_covgs
        self.out_json["phylogenetics"]["species"] = self.species_covgs
        self.out_json["phylogenetics"]["lineage"] = self.lineage_covgs
        if not self.verbose:
            self.out_json["phylogenetics"] = self.choose_best(
                self.out_json["phylogenetics"])
        self._add_unknown_where_empty(
            self.out_json["phylogenetics"]["phylo_group"])
        self._add_unknown_where_empty(
            self.out_json["phylogenetics"]["sub_complex"])
        self._add_unknown_where_empty(
            self.out_json["phylogenetics"]["species"])
        self._add_unknown_where_empty(
            self.out_json["phylogenetics"]["lineage"])

    def _bases_covered(self, percent_coverage, length):
        return sum([percent_coverage[i] * length[i]
                    for i in range(len(length))])

    def _aggregate(self, covgs, threshold=5):
        del_phylo_groups = []
        for phylo_group, covg_dict in covgs.items():
            percent_coverage = covg_dict["percent_coverage"]
            length = covg_dict["length"]
            bases_covered = self._bases_covered(percent_coverage, length)
            total_bases = covg_dict["total_bases"]
            total_percent_covered = round(bases_covered / total_bases, 3)
            _median = covg_dict.get("median", [0])
            minimum_percentage_coverage_required = percent_coverage_from_expected_coverage(
                self.expected_depth) * self.threshold.get(phylo_group, DEFAULT_THRESHOLD)
            if total_percent_covered < minimum_percentage_coverage_required or median(
                    _median) < 0.1 * self.expected_depth:
                # Remove low coverage nodes
                _index = [
                    i for i,
                    d in enumerate(_median) if d > 0.1 *
                    self.expected_depth]
                percent_coverage = [percent_coverage[i] for i in _index]
                length = [length[i] for i in _index]
                bases_covered = self._bases_covered(percent_coverage, length)
                _median = [_median[i] for i in _index]
                total_percent_covered = round(bases_covered / total_bases, 3)
            if total_percent_covered > threshold:
                if phylo_group == "Mycobacterium_llatzerense":  # Mistake in panel
                    phylo_group = "Mycobacterium_mucogenicum"
                covgs[phylo_group] = {
                    "percent_coverage": total_percent_covered,
                    "median_depth": median(_median)}
            else:
                del_phylo_groups.append(phylo_group)
        for phylo_group in del_phylo_groups:
            del covgs[phylo_group]

    def choose_best(self, phylogenetics):
        # Get all the phylo_groups present.
        phylo_groups = self._get_present_phylo_groups(
            phylogenetics["phylo_group"])
        phylogenetics["phylo_group"] = phylo_groups
        sub_complexes = self._get_present_phylo_groups(
            phylogenetics["sub_complex"],
            mix_threshold=90)
        phylogenetics["sub_complex"] = sub_complexes
        # for each phylo_group, get the best species in the phylo_group (using
        # sub_complex info where possible)
        species = {}
        for pg in phylo_groups.keys():
            if self.hierarchy:
                allowed_species = flatten([self.hierarchy.dict[pg]["children"][subc]["children"].keys(
                ) for subc in self.hierarchy.dict[pg]["children"].keys() if subc != "Unknown"])
                species_to_consider = {k: phylogenetics["species"].get(
                    k, {"percent_coverage": 0}) for k in allowed_species}
            else:
                species_to_consider = phylogenetics["species"]
            best_species = self._get_present_phylo_groups(
                species_to_consider,
                mix_threshold=90)
            species.update(best_species)
        phylogenetics["species"] = species
        # For each species, get the best sub species where applicable
        sub_species = {}
        for s in species.keys():
            if self.hierarchy:
                allowed_sub_species = self.hierarchy.get_children(s)
                sub_species_to_consider = {k: phylogenetics["lineage"].get(
                    k, {"percent_coverage": 0}) for k in allowed_sub_species}
            else:
                sub_species_to_consider = phylogenetics.get("lineage", {})
            best_sub_species = self._get_best_coverage_dict(
                sub_species_to_consider)
            sub_species.update(best_sub_species)
        phylogenetics["lineage"] = sub_species
        return phylogenetics

    def _get_present_phylo_groups(self, phylo_groups, mix_threshold=50):
        # If there are more than one pg, if there are more that one above high
        # conf threshold return both
        if not phylo_groups:
            return phylo_groups
        high_confidence_phylo_groups = [
            pg for pg,
            d in phylo_groups.items() if d["percent_coverage"] > mix_threshold]
        if len(high_confidence_phylo_groups) > 1:
            # high_confidence_phylo_groups
            return {k: phylo_groups.get(k, {"percent_coverage": 0})
                    for k in high_confidence_phylo_groups}
        else:
            # Otherwise return best hit
            return self._get_best_coverage_dict(phylo_groups)

    def _get_best_coverage_dict(self, coverage_dict):
        # If there are more than one pg, if there are more that one above high
        # conf threshold return both
        if not coverage_dict:
            return coverage_dict
        sorted_coverage_dict = sorted(
            coverage_dict.items(),
            key=lambda x: (x[1]["percent_coverage"], x[
                           1].get("median_depth", 0)),
            reverse=True)
        if (sorted_coverage_dict[0][1]["percent_coverage"]) > 0:
            return {sorted_coverage_dict[0][0]: sorted_coverage_dict[0][1]}
        else:
            return {}


class AMRSpeciesPredictor(SpeciesPredictor):

    def __init__(
            self,
            phylo_group_covgs,
            sub_complex_covgs,
            species_covgs,
            lineage_covgs,
            verbose=False,
            hierarchy_json_file=None):
        super(
            AMRSpeciesPredictor,
            self).__init__(
            phylo_group_covgs,
            sub_complex_covgs,
            species_covgs,
            lineage_covgs,
            verbose=verbose,
            hierarchy_json_file=hierarchy_json_file)

    def is_saureus_present(self):
        return "Staphaureus" in self.out_json["phylogenetics"]["phylo_group"]

    def is_mtbc_present(self):
        return "Mycobacterium_tuberculosis_complex" in self.out_json[
            "phylogenetics"]["phylo_group"]

    def is_ntm_present(self):
        return "Non_tuberculosis_mycobacterium_complex" in self.out_json[
            "phylogenetics"]["phylo_group"]

    def is_gram_neg_present(self):
        return self.is_klebsiella_pneumoniae_present(
        ) or self.is_escherichia_coli_present()

    def is_klebsiella_pneumoniae_present(self):
        return "Klebsiella_pneumoniae" in self.out_json[
            "phylogenetics"]["species"]

    def is_escherichia_coli_present(self):
        return "Escherichia_coli" in self.out_json["phylogenetics"]["species"]

    def contamination_depths(self):
        contamination_depths = []
        ignore = []
        if self.is_saureus_present():
            ignore.append("Saureus")
        elif self.is_mtbc_present():
            ignore.append("Mtuberculosis")
        elif self.is_escherichia_coli_present():
            ignore.append("Escherichia_coli")
        elif self.is_klebsiella_pneumoniae_present():
            ignore.append("Klebsiella_pneumoniae")
        for node, covg_collection in self.species_covgs.items():
            if node not in ignore:
                contamination_depths.append(covg_collection["median_depth"])
        return contamination_depths
