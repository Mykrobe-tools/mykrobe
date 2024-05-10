import os
import pytest

from mykrobe.probes.models import Mutation
from mykrobe.annotation.genes import GeneAminoAcidChangeToDNAVariants


@pytest.fixture(scope="module")
def ref_data():
    data_dir = os.path.join("tests", "ref_data")
    reference_filepath = f"{data_dir}/NC_000962.3.fasta"
    reference = os.path.basename(reference_filepath).split(".fa")[0]
    aa2dna = GeneAminoAcidChangeToDNAVariants(
        f"{data_dir}/NC_000962.3.fasta", f"{data_dir}/NC_000962.3.gb"
    )
    return reference, aa2dna


def test_mutation_name_forward_strand(ref_data):
    reference, aa2dna = ref_data
    gene = "rpoB"
    mutation_string = "S450L"
    is_protein_coding_var = True
    assert set(
        aa2dna.get_variant_names(gene, mutation_string, is_protein_coding_var)
    ) == set(
        [
            "TCG761154TTA",
            "TCG761154TTG",
            "TCG761154CTA",
            "TCG761154CTT",
            "TCG761154CTC",
            "TCG761154CTG",
        ]
    )
    mutation = Mutation(
        reference=reference,
        var_name="TCG761154TTA",
        gene=aa2dna.get_gene("rpoB"),
        mut="S450L",
    )
    assert mutation.mutation_output_name == "S450L"


def test_mutation_name_reverse_strand(ref_data):
    reference, aa2dna = ref_data
    gene = "gid"
    mutation_string = "I11N"
    is_protein_coding_var = True
    assert set(
        aa2dna.get_variant_names(gene, mutation_string, is_protein_coding_var)
    ) == set(["GAT4408170ATT", "GAT4408170GTT"])
    mutation = Mutation(
        reference=reference,
        var_name="GAT4408170ATT",
        gene=aa2dna.get_gene("gid"),
        mut="I11N",
    )
    assert mutation.mutation_output_name == "I11N"


def test_mutation_name_dna_space(ref_data):
    reference, aa2dna = ref_data
    gene = "pncA"
    mutation_string = "C18CCA"
    is_protein_coding_var = False
    assert set(
        aa2dna.get_variant_names(gene, mutation_string, is_protein_coding_var)
    ) == set(["G2289224TGG"])
    mutation = Mutation(
        reference=reference,
        var_name=aa2dna.get_variant_names(gene, mutation_string, is_protein_coding_var)[
            0
        ],
        gene=aa2dna.get_gene(gene),
        mut=mutation_string,
    )
    assert mutation.mutation_output_name == "C18CCA"
