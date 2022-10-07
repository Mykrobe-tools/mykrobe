#!/usr/bin/env python3
"""This script can be used to combine multiple mykrobe CSV/JSON output files into a
a sheet for easy comparison. The columns are drugs, the rows are samples, and the cells
are the susceptibility prediction.

USAGE:

# combine two CSV reports into a file called samples.csv
combine.py -o samples.csv sample1.csv sample2.csv

# combine a JSON report with a CSV report
combine.py -o samples.csv sample1.json sample2.csv

# combine all JSON and CSV files in a directory
combine.py -o samples.csv reports/
"""
import json
import sys
import argparse
from pathlib import Path
from typing import Dict


def extract_csv_data(path: Path) -> Dict[str, str]:
    report = dict()
    with open(path) as fp:
        header_fields = [s.strip('"') for s in next(fp).split(",")]
        sample_idx = header_fields.index("sample")
        prediction_idx = header_fields.index("susceptibility")
        drug_idx = header_fields.index("drug")
        species_idx = header_fields.index("species")
        lineage_idx = header_fields.index("lineage")

        samples = set()
        for row in fp:
            fields = [s.strip('"') for s in row.split(",")]
            sample = fields[sample_idx]
            samples.add(sample)
            drug = fields[drug_idx]
            pred = fields[prediction_idx]
            report["species"] = fields[species_idx]
            report["lineage"] = fields[lineage_idx]
            report[drug] = pred

        if len(samples) > 1:
            eprint(f"ERROR: Got more than one sample in {path}")
            sys.exit(1)

        report["sample"] = sample
        print(report)
        return report


def extract_json_data(path: Path) -> Dict[str, str]:
    report = dict()
    with open(path) as fp:
        data = json.load(fp)

    samples = list(data.keys())
    if len(samples) > 1:
        eprint(f"ERROR: Got moe than one sample in {path}")
        sys.exit(1)

    report["sample"] = samples[0]

    for drug, info in data[samples[0]]["susceptibility"].items():
        pred = info["predict"]
        report[drug] = pred

    phylo = data[samples[0]].get("phylogenetics", {})
    report["species"] = ";".join(list(phylo.get("species", {}).keys()))
    lineages = [
        s.replace("lineage", "") for s in phylo.get("lineage", {}).get("lineage", {})
    ]
    report["lineage"] = ";".join(lineages)

    return report


def extract_data(path: Path) -> Dict[str, str]:
    if path.suffix == ".csv":
        return extract_csv_data(path)
    elif path.suffix == ".json":
        return extract_json_data(path)
    else:
        eprint(f"ERROR: {path} is an unknown file type")
        sys.exit(1)


def eprint(msg: str):
    print(msg, file=sys.stderr)


def main():
    if sys.version_info < (3, 6, 0):
        eprint("You need python 3.6 or later to run this script")
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "src",
        type=Path,
        metavar="DIR/FILES",
        nargs="+",
        help="A directory containing mykrobe reports, or a list of mykrobe report files",
    )
    parser.add_argument(
        "-d",
        "--delim",
        help="Delimiting character to use in the output file [default: %(default)s]",
        default=",",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="File path to save output to [default: stdout]",
        type=Path,
    )

    args = parser.parse_args()

    delim = args.delim
    input_files = []
    for p in args.src:
        if p.is_dir():
            input_files.append(list(p.glob("*.csv")))
            input_files.append(list(p.glob("*.json")))
        else:
            input_files.append(p)

    eprint(f"Found {len(input_files)} mykrobe reports")

    drugs = set()
    reports = dict()
    for p in input_files:
        data = extract_data(p)
        reports[p] = data
        for k in data:
            if k in ("sample", "species", "lineage"):
                continue
            drugs.add(k)

    fp_out = args.output.open("w") if args.output is not None else sys.stdout
    print(delim.join(["sample", "species", "lineage", *sorted(drugs)]), file=fp_out)

    for p, data in reports.items():
        row = []
        for k in ["sample", "species", "lineage"]:
            row.append(data[k])

        for drug in sorted(drugs):
            row.append(data.get(drug, ""))

        print(delim.join(row), file=fp_out)

    fp_out.close()


if __name__ == "__main__":
    main()
