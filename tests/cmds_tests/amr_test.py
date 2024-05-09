import copy

from mykrobe import amr


def test_add_ncbi_species_names_to_phylo_dict():
    phylo_expect = {
        "phylo_group": {
            "complex_name": {
                "percent_coverage": 99,
                "median_depth": 10,
            },
        },
        "sub_complex": {
            "sub_complex_name": {"foo": "bar"},
        },
    }

    phylo = copy.deepcopy(phylo_expect)
    amr.add_ncbi_species_names_to_phylo_dict(phylo, {})
    assert phylo == phylo_expect

    phylo_expect["species"] = {
        "s1": {
            "percent_coverage": 99,
            "median_depth": 10,
        },
        "s2": {
            "percent_coverage": 98,
            "median_depth": 11,
        },
    }
    phylo = copy.deepcopy(phylo_expect)
    amr.add_ncbi_species_names_to_phylo_dict(phylo, {})
    phylo_expect["species"]["s1"]["ncbi_names"] = "UNKNOWN"
    phylo_expect["species"]["s2"]["ncbi_names"] = "UNKNOWN"
    assert phylo == phylo_expect

    ncbi_names = {"s1": "s1_ncbi_name"}
    amr.add_ncbi_species_names_to_phylo_dict(phylo, ncbi_names)
    phylo_expect["species"]["s1"]["ncbi_names"] = "s1_ncbi_name"
    assert phylo == phylo_expect

    ncbi_names["s2"] = "s2_ncbi_name"
    amr.add_ncbi_species_names_to_phylo_dict(phylo, ncbi_names)
    phylo_expect["species"]["s2"]["ncbi_names"] = "d2_ncbi_name"
