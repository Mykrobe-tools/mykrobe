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
    for k, v in d.get("phylogenetics", {}).get("lineage", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_file_name(f):
    sample_name = os.path.basename(f).split('.')[0]
    return sample_name


def get_sample_name(f):
    return list(f.keys())[0]

def get_plate_name(f):
    return f.split('/')[-3]

def get_expected_depth(d):
    return str(d.get("expected_depth", -1))


def get_mean_read_length(d):
    return str(d.get("mean_read_length", -1))


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
            wt_depth = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("median_depth")
            alt_depth = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("median_depth")

            wt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("percent_coverage")
            alt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("percent_coverage")
            if wt_per_cov < 100:
                wt_depth = 0
            if alt_per_cov <100:
                alt_depth =0 
            conf = variant_call.get('info',{}).get('conf',"1")

            variants.append(":".join([variant_name,
                             str(int(alt_depth)),str(int(wt_depth)), str(int(conf))
                           ]))
    return ";".join(variants)


def json_to_csv(json_dict):
    header = [
        "drug",
        "susceptibility",  
        "variants (dna_variant-AA_variant:alt_depth:ref_depth:conf) [use --format json for more info]",
        "genes (prot_mut-ref_mut:percent_covg:depth) [use --format json for more info]",              
        "sample",

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
            row = [drug,
                call.get(
                    "predict",
                    'N'),
                called_by_variants, 
                called_by_genes,  
                sample_name,
                              
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