import json
import os
import pytest
import subprocess

from mykrobe.parser import parser
import mykrobe


def test_get_panels_and_run_predict_on_test_reads():
    # This is a basic end-to-end test using the Mtb test reads that are on
    # figshare. We'll run mykrobe predict on the reads, and check we get
    # INH resistance, and that the lineage call is 4.10.
    # This means downloading the panel data, downloading the reads, running
    # mykrobe predict, and checking the relevant parts of the JSON output file.
    # Note: not calling this via subprocess because we want to get an
    # accurate coverage report.
    out_dir = "test.myk_tb_end_to_end"
    subprocess.check_output(f"rm -rf {out_dir}", shell=True)
    os.mkdir(out_dir)

    # Test we can download the panels ok, and basic check we got TB and the
    # panel we want to use
    panels_dir = os.path.join(out_dir, "panels")
    args = parser.parse_args(["panels", "update_metadata", "--panels_dir", panels_dir])
    mykrobe.cmds.panels.update_metadata(parser, args)
    args = parser.parse_args(
        ["panels", "update_species", "--panels_dir", panels_dir, "all"]
    )
    mykrobe.cmds.panels.update_species(parser, args)
    panels_json = os.path.join(panels_dir, "manifest.json")
    assert os.path.exists(panels_json)
    with open(panels_json) as f:
        panels_manifest = json.load(f)
    assert "tb" in panels_manifest
    tb_json = os.path.join(panels_dir, "tb", "manifest.json")
    assert os.path.exists(tb_json)
    with open(tb_json) as f:
        tb_manifest = json.load(f)
    panel_to_use = "202010"
    assert "panels" in tb_manifest
    assert panel_to_use in tb_manifest["panels"]

    # Get the test reads from figshare and run predict on them
    reads_fq = os.path.join(out_dir, "reads.fq.gz")
    predict_json = os.path.join(out_dir, "predict.json")
    if os.path.exists(predict_json):
        os.unlink(predict_json)
    subprocess.check_output(
        f"wget -q -O {reads_fq} https://ndownloader.figshare.com/files/21059229",
        shell=True,
    )
    sample_name = "TEST_SAMPLE"
    # Important! We need to use a different name for the skeletons directory,
    # because the default mykrobe/data/skeletons breaks with the coverage
    # module, resulting in no coverage report.
    skeleton_dir = os.path.join(out_dir, "skeletons")
    args = parser.parse_args(
        [
            "predict",
            "--panel",
            panel_to_use,
            "--skeleton_dir",
            skeleton_dir,
            "--sample",
            sample_name,
            "--species",
            "tb",
            "--format",
            "json",
            "--output",
            predict_json,
            "--seq",
            reads_fq,
            "--panels_dir",
            panels_dir,
        ]
    )
    mykrobe.cmds.amr.run(parser, args)

    # Basic checks of the results:
    # - we got a JSON file
    # - there's one sample in the json, with the expected name
    # - all drugs susceptible, except isoniazid should be resistant
    # - lineage is 4.10.
    with open(predict_json) as f:
        amr_results = json.load(f)
    assert len(amr_results) == 1
    assert sample_name in amr_results
    for drug, results in amr_results[sample_name]["susceptibility"].items():
        if drug == "Isoniazid":
            assert results["predict"] == "R"
        else:
            assert results["predict"] == "S"
    assert amr_results[sample_name]["phylogenetics"]["lineage"]["lineage"] == [
        "lineage4.10"
    ]
    subprocess.check_output(f"rm -rf {out_dir}", shell=True)
