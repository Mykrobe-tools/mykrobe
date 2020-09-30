import os
import pytest
from mykrobe.probes.models import Mutation
from mykrobe.probes import load_dna_vars_txt_file

def test_load_dna_vars_txt_file():
    infile = os.path.join("tests", "probe_tests", "test_probe_generation.load_dna_vars_txt_file.tsv")
    got_mutations, got_lineage = load_dna_vars_txt_file(infile, "ref")
    expect_mutations = [
        Mutation(reference="ref", var_name="G42A"),
        Mutation(reference="ref", var_name="C52G"),
        Mutation(reference="ref", var_name="C62G"),
        Mutation(reference="ref", var_name="A72T"),
    ]
    expect_lineage = {
            "G42A": {"name": "lineage1", "use_ref_allele": False},
            "C52G": {"name": "lineage1.2", "use_ref_allele": True},
            "A72T": {"name": "lineage1.2.3", "use_ref_allele": False, "report_name": "lineage1.2.3_report_name"},
    }
    assert got_mutations == expect_mutations
    assert got_lineage == expect_lineage

