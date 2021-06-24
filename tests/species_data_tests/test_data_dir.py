import copy
import json
import os
import pytest
import shutil
import time

from mykrobe.species_data import DataDir

data_dir = os.path.dirname(os.path.abspath(__file__))

def test_data_dir():
    # This is a long test that runs through the entire reference data process
    # that a user might do: update metadata, install a species, update a
    # species, and remove a species. Checks along the way that all the metadata
    # etc is correct.

    # Create an empty data dir
    temp_dir = "tmp.test_data_dir"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    ddir = DataDir(temp_dir)
    assert ddir.manifest == {}
    ddir.create_root()
    assert os.path.exists(temp_dir)
    assert not ddir.is_locked()
    ddir.start_lock()
    assert ddir.is_locked()
    ddir.stop_lock()
    assert not ddir.is_locked()

    # Update the manifest, which has species available for installation
    species1_tarball = os.path.join(data_dir, "species1_data.20200101.tar.gz")
    manifest_data = {
        "species1": {"version": "20200101", "url": species1_tarball},
        "species2": {"version": "20190211", "url": "species2_url"},
    }
    manifest_json = "tmp.species_data_test.json"
    with open(manifest_json, "w") as f:
        json.dump(manifest_data, f)
    ddir.update_manifest(filename=manifest_json)
    os.unlink(manifest_json)
    expect_manifest = {
        "species1": {"installed": None, "latest": copy.copy(manifest_data["species1"])},
        "species2": {"installed": None, "latest": copy.copy(manifest_data["species2"])},
    }
    assert ddir.manifest == expect_manifest
    assert not ddir.is_locked()
    assert ddir.all_species_list() == ["species1", "species2"]
    assert ddir.installed_species() == []
    assert not ddir.species_is_installed("species1")
    assert not ddir.species_is_installed("species2")
    assert ddir.get_species_dir("species1") is None
    assert ddir.get_species_dir("species2") is None
    with pytest.raises(ValueError):
        ddir.get_species_dir("unknown species")

    # Add species1 from a tarball file on disk, and check everything looks
    # correctly updated.
    ddir.add_or_replace_species_data(species1_tarball)
    assert not ddir.is_locked()
    assert ddir.species_is_installed("species1")
    assert ddir.installed_species() == ["species1"]
    expect_manifest_with_species1 = copy.deepcopy(expect_manifest)
    expect_manifest_with_species1["species1"]["installed"] = copy.copy(manifest_data["species1"])
    assert ddir.manifest == expect_manifest_with_species1

    # Want to test we can get species1 from the data dir.
    # But getting filesystem delay issues with the species1 manifest.json file,
    # which is loaded by ddir.get_species_dir("species1") below. Python
    # thinks it's there according to os.path.exists(), but then opening the
    # file throws a FileNotFoundError?! Hence the next loop.
    for i in range(3):
        time.sleep(0.5)
        try:
            species1 = ddir.get_species_dir("species1")
        except:
            pass
        if species1 is not None:
            break
    assert species1 is not None

    # Make a new metadata manifest, which has a newer version of species1,
    # and use it to update the data directory.
    species1_tarball = os.path.join(data_dir, "species1_data.20200801.tar.gz")
    manifest_data = {
        "species1": {"version": "20200801", "url": species1_tarball},
        "species2": {"version": "20190211", "url": "species2_url"},
    }
    manifest_json = "tmp.species_data_test.json"
    with open(manifest_json, "w") as f:
        json.dump(manifest_data, f)
    ddir.update_manifest(filename=manifest_json)
    os.unlink(manifest_json)

    # If we try to update species1 again, should fail because it's already
    # installed and we didn't force it. But should not have removed the existing
    # install of species1.
    with pytest.raises(RuntimeError):
        ddir.add_or_replace_species_data(species1_tarball)
    assert ddir.is_locked()
    ddir.stop_lock()
    assert ddir.species_is_installed("species1")

    # Now update species1 with the force option and check it worked.
    ddir.add_or_replace_species_data(species1_tarball, force=True)
    expect_manifest_with_species1["species1"]["latest"] = copy.copy(manifest_data["species1"])
    expect_manifest_with_species1["species1"]["installed"] = copy.copy(manifest_data["species1"])
    assert ddir.manifest == expect_manifest_with_species1
    assert ddir.species_is_installed("species1")
    assert ddir.installed_species() == ["species1"]

    # Test removing species. Unknown species should fail, species1 should work
    with pytest.raises(ValueError):
        ddir.remove_species("unknown species")
    ddir.remove_species("species1")
    assert not ddir.species_is_installed("species1")
    assert ddir.installed_species() == []
    expect_manifest["species1"]["latest"] = copy.copy(manifest_data["species1"])
    assert ddir.manifest == expect_manifest
    shutil.rmtree(temp_dir)
