import copy
import json
import pytest

import anytree
from anytree.exporter import JsonExporter

from mykrobe.metagenomics import LineagePredictor


def test_lineage_to_itself_plus_parents():
    f = LineagePredictor._lineage_to_itself_plus_parents
    assert f("a") == ["a"]
    assert f("a.1") == ["a", "a.1"]
    assert f("a.1.2") == ["a", "a.1", "a.1.2"]

def test_constructor_makes_tree():
    # Note: deliberately miss out lineage1.1 to check that node
    # still gets made. We're allowing the user to not specify every single
    # node in their lineage scheme
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": False},
        "var1a": {"name": "lineage1", "use_ref_allele": False},
        "var1.2": {"name": "lineage1.2", "use_ref_allele": False},
        "var1.1.1": {"name": "lineage1.1.1", "use_ref_allele": False},
        "var2": {"name": "lineage2", "use_ref_allele": False},
        "var2.1": {"name": "lineage2.1", "use_ref_allele": False},
    }
    lin_pred = LineagePredictor(variant_to_lineage)
    exporter = JsonExporter()
    got_tree = json.loads(exporter.export(lin_pred.tree_root))
    expect_tree = {
        "name": "root",
        "children": [
            {
                "name": "lineage1",
                "children": [
                    {"name": "lineage1.2"},
                    {"name": "lineage1.1", "children": [{"name": "lineage1.1.1"}],},
                ],
            },
            {"name": "lineage2", "children": [{"name": "lineage2.1"}],},
        ],
    }
    assert got_tree == expect_tree


def test_score_each_lineage_node():
    lineage_calls = {
        "lineage1": {
            "var1": {"genotype": [0, 0], "info": {"filter": [], "conf": 500},},
            "var1a": {"genotype": [1, 1], "info": {"filter": [], "conf": 1000},},
        },
        "lineage1.1": {
            "var1.1": {"genotype": [1, 1], "info": {"filter": [], "conf": 20},},
        },
    }
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": True},
        "var1a": {"name": "lineage1", "use_ref_allele": True},
        "var1.1": {"name": "lineage1.1", "use_ref_allele": False},
    }
    expect = {"lineage1": 500, "lineage1.1": 20}
    lin_pred = LineagePredictor(variant_to_lineage)
    assert lin_pred._score_each_lineage_node(lineage_calls) == expect


def test_get_paths_and_scores():
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": True},
        "var1a": {"name": "lineage1", "use_ref_allele": True},
        "var1.1": {"name": "lineage1.1", "use_ref_allele": False},
        "var1.2": {"name": "lineage1.2", "use_ref_allele": False},
        "var1.1.1": {"name": "lineage1.1.1", "use_ref_allele": False},
        "var2": {"name": "lineage2", "use_ref_allele": False},
        "var2.1": {"name": "lineage2.1", "use_ref_allele": False},
        "var2.2": {"name": "lineage2.2", "use_ref_allele": False},
    }

    lineage_calls = {
        "lineage1": {
            "var1": {"genotype": [0, 0], "info": {"filter": [], "conf": 500},},
            "var1a": {"genotype": [1, 1], "info": {"filter": [], "conf": 1000},},
        },
        "lineage1.1": {
            "var1.1": {"genotype": [1, 1], "info": {"filter": [], "conf": 1000},},
        },
        "lineage1.1.1": {
            "var1.1.1": {"genotype": [1, 1], "info": {"filter": [], "conf": 100},},
        },
        "lineage2": {"var2": {"genotype": [1, 1], "info": {"filter": [], "conf": 1},},},
        "lineage2.1": {
            "var2": {"genotype": [0, 0], "info": {"filter": [], "conf": 100},},
        },
    }

    lin_pred = LineagePredictor(variant_to_lineage)
    got = lin_pred._get_paths_and_scores(lineage_calls)
    expect = [
        {
            "score": 1600,
            "lineage": "lineage1.1.1",
            "scores": {"lineage1": 500, "lineage1.1": 1000, "lineage1.1.1": 100},
        },
        {"score": 500, "lineage": "lineage1", "scores": {"lineage1": 500}},
        {"score": 1, "lineage": "lineage2", "scores": {"lineage2": 1}},
    ]
    assert got == expect


def test_call_lineage_using_conf_scores():
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": False},
        "var1a": {"name": "lineage1", "use_ref_allele": False},
        "var1.1": {"name": "lineage1.1", "use_ref_allele": False},
        "var1.1.1": {"name": "lineage1.1.1", "use_ref_allele": False},
        "var2": {"name": "lineage2", "use_ref_allele": False},
        "var2.1": {"name": "lineage2.1", "use_ref_allele": False},
    }
    lin_pred = LineagePredictor(variant_to_lineage)
    assert lin_pred.call_lineage_using_conf_scores({}) is None

    var1_call = {"genotype": [0, 0], "info": {"filter": [], "conf": 200}}
    var1a_call = {"genotype": [1, 1], "info": {"filter": [], "conf": 300}}
    var1_1_call = {"genotype": [1, 1], "info": {"filter": [], "conf": 1000}}
    var1_1_1_call = {"genotype": [1, 1], "info": {"filter": [], "conf": 100}}
    var2_call = {"genotype": [1, 1], "info": {"filter": [], "conf": 900}}
    var2_1_call = {"genotype": [1, 1], "info": {"filter": [], "conf": 500}}

    lineage_calls = {
        "lineage1": {"var1": var1_call, "var1a": var1a_call,},
        "lineage1.1": {"var1.1": var1_1_call,},
        "lineage1.1.1": {"var1.1.1": var1_1_1_call,},
    }

    expect = {
        "lineage": ["lineage1.1.1"],
        "calls": {"lineage1.1.1": copy.deepcopy(lineage_calls)},
    }
    assert lin_pred.call_lineage_using_conf_scores(lineage_calls) == expect

    lineage_calls["lineage2"] = {"var2": var2_call}
    lineage_calls["lineage2.1"] = {"var2.1": var2_1_call}
    expect["lineage"].append("lineage2.1")
    expect["calls"]["lineage2.1"] = {
        "lineage2": {"var2": var2_call},
        "lineage2.1": {"var2.1": var2_1_call},
    }
    assert lin_pred.call_lineage_using_conf_scores(lineage_calls) == expect


def test_genotype_each_lineage_node():
    lineage_calls = {
        "lineage1": {"var1": {"genotype": [0, 0]}, "var1a": {"genotype": [1, 1]},},
        "lineage1.1": {"var1.1": {"genotype": [0, 1]},},
    }
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": True},
        "var1a": {"name": "lineage1", "use_ref_allele": True},
        "var1.1": {"name": "lineage1.1", "use_ref_allele": False},
    }
    expect = {"lineage1": 1, "lineage1.1": 0.5}
    lin_pred = LineagePredictor(variant_to_lineage)
    assert lin_pred._genotype_each_lineage_node(lineage_calls) == expect


def test_get_good_paths_using_genotype_calls():
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": True},
        "var1a": {"name": "lineage1", "use_ref_allele": True},
        "var1.1": {"name": "lineage1.1", "use_ref_allele": False},
        "var1.2": {"name": "lineage1.2", "use_ref_allele": False},
        "var1.1.1": {"name": "lineage1.1.1", "use_ref_allele": False},
        "var2": {"name": "lineage2", "use_ref_allele": False},
        "var2.1": {"name": "lineage2.1", "use_ref_allele": False},
        "var2.2": {"name": "lineage2.2", "use_ref_allele": False},
    }

    lineage_calls = {
        "lineage1": {"var1": {"genotype": [0, 0],}, "var1a": {"genotype": [1, 1]},},
        "lineage1.1": {"var1.1": {"genotype": [1, 1]},},
        "lineage1.1.1": {"var1.1.1": {"genotype": [1, 1]},},
        "lineage2": {"var2": {"genotype": [1, 1]},},
        "lineage2.1": {"var2": {"genotype": [0, 0]},},
    }

    lin_pred = LineagePredictor(variant_to_lineage)
    got = lin_pred._get_good_paths_using_genotype_calls(lineage_calls)
    expect = {
        "lineage1.1.1":
        {
            "genotypes": {"lineage1": 1, "lineage1.1": 1, "lineage1.1.1": 1},
            "good_nodes": 3,
            "tree_depth": 3,
        },
        "lineage2":
        {
            "genotypes": {"lineage2": 1},
            "good_nodes": 1,
            "tree_depth": 1,
        },
    }
    assert got == expect


def test_call_lineage():
    variant_to_lineage = {
        "var1": {"name": "lineage1", "use_ref_allele": False},
        "var1a": {"name": "lineage1", "use_ref_allele": False},
        "var1.1": {"name": "lineage1.1", "use_ref_allele": False},
        "var1.1.1": {"name": "lineage1.1.1", "use_ref_allele": False, "report_name": "l1-1-1"},
        "var2": {"name": "lineage2", "use_ref_allele": False},
        "var2.1": {"name": "lineage2.1", "use_ref_allele": False},
    }
    lin_pred = LineagePredictor(variant_to_lineage)
    assert lin_pred.call_lineage({}) is None

    var1_call = {"genotype": [0, 0]}
    var1a_call = {"genotype": [1, 1]}
    var1_1_call = {"genotype": [1, 1]}
    var1_1_1_call = {"genotype": [1, 1]}
    var2_call = {"genotype": [1, 1]}
    var2_1_call = {"genotype": [1, 1]}

    lineage_calls = {
        "lineage1": {"var1": var1_call, "var1a": var1a_call,},
        "lineage1.1": {"var1.1": var1_1_call,},
        "lineage1.1.1": {"var1.1.1": var1_1_1_call,},
    }

    expect_lineage_calls = copy.deepcopy(lineage_calls)
    expect_lineage_calls["l1-1-1"] = {"var1.1.1": var1_1_1_call}
    del expect_lineage_calls["lineage1.1.1"]

    expect = {
        "lineage": ["l1-1-1"],
        "calls": {"l1-1-1": expect_lineage_calls},
        "calls_summary": {
            "l1-1-1": {
                "genotypes": {"lineage1": 1, "lineage1.1": 1, "l1-1-1": 1},
                "good_nodes": 3,
                "tree_depth": 3,
            },
        }
    }
    assert lin_pred.call_lineage(lineage_calls) == expect

    lineage_calls["lineage2"] = {"var2": var2_call}
    lineage_calls["lineage2.1"] = {"var2.1": var2_1_call}
    expect["lineage"].append("lineage2.1")
    expect["calls"]["lineage2.1"] = {
        "lineage2": {"var2": var2_call},
        "lineage2.1": {"var2.1": var2_1_call},
    }
    expect["calls_summary"]["lineage2.1"] = {
        "genotypes": {'lineage2': 1, 'lineage2.1': 1},
        "good_nodes": 2,
        "tree_depth": 2,
    }
    assert lin_pred.call_lineage(lineage_calls) == expect
