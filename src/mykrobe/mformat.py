import io
import json
import csv
import os

def get_drugs(drug_list):
    drugs = []
    for f in args.files:
        try:
            d = load_json(f)
        except ValueError:
            d = {}
        for drug in drug_list:
            if drug not in drugs:
                drugs.append(drug)
    return drugs


def get_phylo_group_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("phylo_group", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_species_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("species", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_lineage_string(d):
    s = []
    depth=[]
    per_cov=[]
    lineage_dict = d.get("phylogenetics", {}).get("lineage", {})
    if len(lineage_dict) == 0:
        return "", "", ""

    if set(lineage_dict.keys()) == {"lineage", "calls_summary", "calls"}:
        return ";".join(lineage_dict.get("lineage", [])), "NA", "NA"

    for k, v in lineage_dict.items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)

def get_called_genes(d, drug=None):
    variants = []
    for variant_name, variant_call in d.items():
        if variant_call.get("_cls") == "Call.SequenceCall":
            per_cov = variant_call.get('info',{}).get('coverage',{}).get("percent_coverage")
            depth = variant_call.get('info',{}).get('coverage',{}).get("median_depth")
            variants.append(":".join([variant_name,
                             str(int(per_cov)),str(int(depth))]))
    return ";".join(variants)


def get_variant_calls(d):
    variants = []
    for variant_name, variant_call in d.items():
        if variant_call.get("_cls") != "Call.SequenceCall":
            ref_kmer_count = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("kmer_count")
            alt_kmer_count = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("kmer_count")
            conf = variant_call.get('info',{}).get('conf',"1")

            variants.append(":".join([variant_name,
                             str(int(ref_kmer_count)),str(int(alt_kmer_count)), str(int(conf))
                           ]))
    return ";".join(variants)


def json_to_csv(json_dict):
    header = [
        "sample",
        "drug",
        "susceptibility",
        "variants (dna_variant-AA_variant:ref_kmer_count:alt_kmer_count:conf) [use --format json for more info]",
        "genes (prot_mut-ref_mut:percent_covg:depth) [use --format json for more info]",

        "mykrobe_version",
        "files",
        "probe_sets",
        "genotype_model",
        "kmer_size",
        "phylo_group",
        "species",
        "lineage",
        "phylo_group_per_covg",
        "species_per_covg",
        "lineage_per_covg",
        "phylo_group_depth",
        "species_depth",
        "lineage_depth"
]
    rows = [header]

    for sample_name, d in json_dict.items():
        phylo_group,phylo_group_per_covg,phylo_group_depth  = get_phylo_group_string(d)
        species,species_per_covg,species_depth  = get_species_string(d)
        lineage,lineage_per_covg,lineage_depth  = get_lineage_string(d)
        files = ";".join(d["files"])
        probe_sets = ";".join(d["probe_sets"])
        genotype_model = d["genotype_model"]
        kmer_size = str(d["kmer"])
        drug_list = sorted(d.get('susceptibility', {}).keys())
        drugs = sorted(drug_list)

        if not drugs:
            drugs = ["NA"]


        for drug in drugs:
            call = d.get('susceptibility', {}).get(drug, {})
            called_by_variants = get_variant_calls(call.get("called_by",{}))
            called_by_genes = get_called_genes(call.get("called_by",{}))
            row = [sample_name,
                drug,
                call.get(
                    "predict",
                    'N'),
                called_by_variants,
                called_by_genes,

                d.get("version",{}).get("mykrobe-predictor","-1"),
                files,
                probe_sets,
                genotype_model,
                kmer_size,
                phylo_group,
                species,
                lineage,
                phylo_group_per_covg,
                species_per_covg,
                lineage_per_covg,
                phylo_group_depth,
                species_depth,
                lineage_depth
                ]
            rows.append(row)
    output = io.StringIO()
    writer = csv.writer(output, quoting=csv.QUOTE_NONNUMERIC)
    for row in rows:
        writer.writerow(row)
    csv_string=output.getvalue()
    return csv_string
