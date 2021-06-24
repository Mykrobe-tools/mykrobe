import json
import os
import pytest
import shutil

from mykrobe.species_data import SpeciesDir


def test_species_dir():
    temp_dir = "tmp.test_species_dir"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    with pytest.raises(FileNotFoundError):
        sdir = SpeciesDir(temp_dir)
    os.mkdir(temp_dir)
    with pytest.raises(FileNotFoundError):
        sdir = SpeciesDir(temp_dir)

    manifest_data = {
        "species_name": "species_name",
        "version": "20200821",
        "default_panel": "panel1",
        "panels": {
            "panel1": {
                "description": "description of panel1",
                "reference_genome": "NC42",
                "species_phylo_group": "species_xyz",
                "fasta_files": ["probes.fa"],
                "kmer": 21,
                "json_files": {
                    "amr": "amr.json",
                    "lineage": None,
                    "hierarchy": None,
                },
            },
        }
    }
    with open(os.path.join(temp_dir, "manifest.json"), "w") as f:
        json.dump(manifest_data, f, indent=2, sort_keys=True)
    sdir = SpeciesDir(temp_dir)
    assert sdir.manifest == manifest_data
    assert sdir.species_name() == manifest_data["species_name"]
    assert sdir.version() == manifest_data["version"]
    assert sdir.panel_name == manifest_data["default_panel"]
    assert sdir.panel_names() == ["panel1"]
    assert sdir.kmer() == 21
    assert sdir.species_phylo_group() == "species_xyz"
    panel1 = manifest_data["panels"]["panel1"]
    assert sdir.reference_genome() == panel1["reference_genome"]
    assert sdir.description() == panel1["description"]
    fasta_files = sdir.fasta_files()
    assert len(fasta_files) == 1
    assert fasta_files == [os.path.abspath(os.path.join(temp_dir, "probes.fa"))]
    with open(fasta_files[0], "w"):
        pass
    amr_json = sdir.json_file("amr")
    assert amr_json == os.path.abspath(os.path.join(temp_dir, "amr.json"))
    with open(amr_json, "w") as f:
        pass
    assert sdir.sanity_check()

    sdir.panel["fasta_files"] = []
    assert not sdir.sanity_check()
    sdir.panel["fasta_files"] = None
    assert not sdir.sanity_check()
    sdir.panel["fasta_files"] = ["probes.oops_typo.fa"]
    assert not sdir.sanity_check()
    sdir.panel["fasta_files"] = ["probes.fa"]
    assert sdir.sanity_check()
    os.unlink(fasta_files[0])
    assert not sdir.sanity_check()
    with open(fasta_files[0], "w"):
        pass
    assert sdir.sanity_check()

    sdir.panel["json_files"]["amr"] = None
    assert not sdir.sanity_check()
    sdir.panel["json_files"]["amr"] = "does not exist.json"
    assert not sdir.sanity_check()
    sdir.panel["json_files"]["amr"] = "amr.json"
    assert sdir.sanity_check()
    os.unlink(amr_json)
    assert not sdir.sanity_check()
    with open(amr_json, "w") as f:
        pass
    assert sdir.sanity_check()


    manifest_data["panels"]["panel2"] = {
        "description": "description of panel2",
        "reference_genome": "NC42",
        "kmer": 31,
        "species_phylo_group": "species_xyz",
        "fasta_files": ["probes2.fa"],
        "json_files": {
            "amr": "amr2.json",
            "lineage": None,
            "hierarchy": None,
        },
    }
    manifest_data["default_panel"] = "panel2"
    with open(os.path.join(temp_dir, "manifest.json"), "w") as f:
        json.dump(manifest_data, f, indent=2, sort_keys=True)
    sdir = SpeciesDir(temp_dir)
    assert sdir.manifest == manifest_data
    assert sdir.panel_name == "panel2"
    assert not sdir.sanity_check()
    fasta_files = sdir.fasta_files()
    assert len(fasta_files) == 1
    assert fasta_files == [os.path.abspath(os.path.join(temp_dir, "probes2.fa"))]
    with open(fasta_files[0], "w"):
        pass
    amr_json = sdir.json_file("amr")
    assert amr_json == os.path.abspath(os.path.join(temp_dir, "amr2.json"))
    with open(amr_json, "w") as f:
        pass
    assert sdir.sanity_check()
    shutil.rmtree(temp_dir)
