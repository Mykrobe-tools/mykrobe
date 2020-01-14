import os
import json
import logging

from mykrobe.utils import unique
from mykrobe.utils import flatten
from mykrobe.utils import get_params

from mykrobe.predict.models import MykrobePredictorSusceptibilityResult

from mykrobe.variants.schema.models import VariantCall
from mykrobe.variants.schema.models import SequenceCall

from pprint import pprint


logger = logging.getLogger(__name__)

DEFAULT_MIN_GENE_CN = 0.03
DEFAULT_MIN_VARIANT_CN = 0.1


def copy_number(call):
    coverage = call.get('info', {}).get("coverage")
    try:
        alternate_depth = coverage.get("alternate").get("median_depth")
        wt_depth = coverage.get("reference").get("median_depth")
    except:
        alternate_depth = coverage.get("median_depth")
        wt_depth = call.get('info', {}).get("expected_depths")[0]

    return round(float(alternate_depth) / (alternate_depth + wt_depth), 2)


def depth_on_alternate(call):
    coverage = call.get('info', {}).get("coverage")
    try:
        alternate_depth = coverage.get("alternate").get("median_depth")
    except:
        alternate_depth = coverage.get("median_depth")
    return alternate_depth


def is_filtered(call):
    info = call.get('info', {})
    _filter = info.get('filter', [])
    if _filter:
        return True
    else:
        return False


class BasePredictor(object):

    def __init__(self, variant_calls, called_genes,
                 base_json={}, depth_threshold=3, ignore_filtered=True,
                 ignore_minor_calls=False,
                 variant_to_resistance_json_fp=None):
        if variant_to_resistance_json_fp:
            self.variant_or_gene_name_to_resistance_drug = load_json(
                variant_to_resistance_json_fp)
        else:
            self.variant_or_gene_name_to_resistance_drug = load_json(
                self.default_variant_to_resistance_drug)
        self.variant_calls = variant_calls
        self.called_genes = called_genes
        self.drugs = self._get_drug_list_from_variant_to_resistance_drug()
        self.resistance_prediction = self._create_initial_resistance_prediction()
        self.out_json = base_json
        self._cn_threshold = {
            "ermA": 0.19,
            "ermB": 0.19,
            "ermC": 0.19,
            "ermT": 0.19,
            "ermY": 0.19,
            "fusA": 0.03,
            "fusC": 0.03,
            "aacAaphD": 0.04,
            "mecA": 0.06,
            "mupA": 0.21,
            "blaZ": 0.04,
            "tetK": 0.13
        }
        self.depth_threshold = depth_threshold
        self.ignore_filtered = ignore_filtered
        self.ignore_minor_calls = ignore_minor_calls

    def _create_initial_resistance_prediction(self):
        self.result = MykrobePredictorSusceptibilityResult(susceptibility = dict(
            (k, {"predict": "N"}) for k in self.drugs))
        self.resistance_predictions = self.result.susceptibility

    def _get_drug_list_from_variant_to_resistance_drug(self):
        return unique(flatten(self.variant_or_gene_name_to_resistance_drug.values()))

    def predict_antibiogram(self):
        for allele_name, variant_call in self.variant_calls.items():
            self._update_resistance_prediction(allele_name, variant_call)
        for name, gene in self.called_genes.items():
            if isinstance(gene, list):
                if len(gene) > 1:
                    logger.warning(
                        "Ambigious gene call. Continuing regardless. ")
                gene = gene[0]
            self._update_resistance_prediction(name, gene)



    def _update_resistance_prediction(self, allele_name, variant_or_gene):
        variant_or_gene_names = self._get_names(allele_name)
        if "variant" in variant_or_gene:
            variant_or_gene_names += self._get_names(variant_or_gene["variant"])
        for name in variant_or_gene_names:
            drugs = self._get_drugs(name)
            resistance_prediction = self._resistance_prediction(
                variant_or_gene, variant_or_gene_names)

            for drug in drugs:
                current_resistance_prediction = self.resistance_predictions[
                    drug]["predict"]
                assert resistance_prediction is not None
                if current_resistance_prediction == "N":
                    self.resistance_predictions[drug][
                        "predict"] = resistance_prediction
                elif current_resistance_prediction in ["I", "S"]:
                    if resistance_prediction in ["r", "R"]:
                        self.resistance_predictions[drug][
                            "predict"] = resistance_prediction
                elif current_resistance_prediction == "r":
                    if resistance_prediction == "R":
                        self.resistance_predictions[drug][
                            "predict"] = resistance_prediction
                if resistance_prediction in ["r", "R"]:
                    variant_or_gene['variant'] = None
                    try:
                        self.resistance_predictions[drug]["called_by"][allele_name] = variant_or_gene
                    except KeyError:
                        self.resistance_predictions[drug]["called_by"] = {}
                        self.resistance_predictions[drug]["called_by"][allele_name] = variant_or_gene

    def _get_names(self, allele_name):
        names = []
        params = get_params(allele_name)
        if params.get("mut"):
            names.append("_".join([params.get("gene"), params.get("mut")]))
        allele_name_split = allele_name.split('?')[0].split('-')
        if len(allele_name_split) > 1:
            names.append(allele_name_split[1])
        else:
            names.append(allele_name_split[0])

        return names

    def _get_drugs(self, name, lower=False):
        if lower:
            name = name.lower()
        try:
            drugs = self.variant_or_gene_name_to_resistance_drug[name]
        except KeyError:
            try:
                drugs = self.variant_or_gene_name_to_resistance_drug[
                    name.split("-")[0]]
            except KeyError:
                talt_name = list(name)
                talt_name[-1] = "X"
                try:
                    drugs = self.variant_or_gene_name_to_resistance_drug[
                        "".join(talt_name)]
                except KeyError:
                    drugs = []
                    if not lower:
                        return self._get_drugs(name, lower=True)
                    else:
                        pass
                        # print ("Warning:NoEntry for %s" % name)
        assert drugs is not None
        return drugs

    def _resistance_prediction(self, variant_or_gene, names):
        __resistance_prediction = None
        if sum(variant_or_gene.get('genotype')) == 2:
            if (is_filtered(variant_or_gene) and self.ignore_filtered) or depth_on_alternate(variant_or_gene) < self.depth_threshold:
                __resistance_prediction = "N"
            elif self._coverage_greater_than_threshold(variant_or_gene, names):
                __resistance_prediction = "R"
            else:
                __resistance_prediction = "S"
        elif sum(variant_or_gene.get('genotype')) == 1:
            if (is_filtered(variant_or_gene) and self.ignore_filtered) or depth_on_alternate(variant_or_gene) < self.depth_threshold:
                __resistance_prediction = "N"
            elif self._coverage_greater_than_threshold(variant_or_gene, names) and not self.ignore_minor_calls:
                __resistance_prediction = "r"
            else:
                __resistance_prediction = "S"
        elif sum(variant_or_gene.get('genotype')) == 0:
            __resistance_prediction = "S"
        else:
            __resistance_prediction = "N"
        return __resistance_prediction

    def _coverage_greater_than_threshold(self, variant_or_gene, names):
        coveage_threshold = DEFAULT_MIN_VARIANT_CN
        for name in names:
            if name in self._cn_threshold:
                coveage_threshold = self._cn_threshold.get(
                    name, DEFAULT_MIN_VARIANT_CN)
        CN_PASS = copy_number(variant_or_gene) > coveage_threshold
        return CN_PASS

    def run(self):
        self.predict_antibiogram()
        return self.result


def load_json(f):
    with open(f, 'r') as infile:
        return json.load(infile)


class TBPredictor(BasePredictor):

    @property
    def default_variant_to_resistance_drug(self):
        self.data_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '../data/predict/tb/'))
        return os.path.join(
            self.data_dir,
            "variant_to_resistance_drug.json")

    def __init__(self, variant_calls, called_genes, base_json={},
                 depth_threshold=3, ignore_filtered=True, ignore_minor_calls=False,
                 variant_to_resistance_json_fp=None):
        super(
            TBPredictor,
            self).__init__(
            variant_calls,
            called_genes,
            base_json,
            depth_threshold=depth_threshold,
            ignore_filtered=ignore_filtered,
                ignore_minor_calls=ignore_minor_calls,
                variant_to_resistance_json_fp=variant_to_resistance_json_fp)


class StaphPredictor(BasePredictor):

    @property
    def default_variant_to_resistance_drug(self):
        self.data_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '../data/predict/staph/'))
        return os.path.join(
            self.data_dir,
            "variant_to_resistance_drug.json")

    def __init__(self, variant_calls, called_genes, base_json={},
                 depth_threshold=3, ignore_filtered=True, ignore_minor_calls=False,
                 variant_to_resistance_json_fp=None):

        super(
            StaphPredictor,
            self).__init__(
            variant_calls,
            called_genes,
            base_json,
            depth_threshold=depth_threshold,
            ignore_filtered=ignore_filtered,
                ignore_minor_calls=ignore_minor_calls,
                variant_to_resistance_json_fp=variant_to_resistance_json_fp)
