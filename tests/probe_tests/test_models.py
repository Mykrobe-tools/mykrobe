import sys
import json
import os

sys.path.append(".")
from mykrobe.probes.models import Mutation
from mykrobe.annotation.genes import GeneAminoAcidChangeToDNAVariants
DATA_DIR = os.path.join("tests", "ref_data")


class TestMutation():

    def setup(self):
        self.reference_filepath = f"{DATA_DIR}/NC_000962.3.fasta"
        self.reference = os.path.basename(
            self.reference_filepath).split('.fa')[0]
        self.aa2dna = GeneAminoAcidChangeToDNAVariants(
            f"{DATA_DIR}/NC_000962.3.fasta",
            f"{DATA_DIR}/NC_000962.3.gb")

    def teardown(self):
        pass

    def test_mutation_name_forward_strand(self):
        gene = "rpoB"
        mutation_string = "S450L"
        is_protein_coding_var = True
        assert set(self.aa2dna.get_variant_names(
            gene, mutation_string, is_protein_coding_var)) == set(["TCG761154TTA", "TCG761154TTG", "TCG761154CTA", "TCG761154CTT", "TCG761154CTC", "TCG761154CTG"])
        mutation = Mutation(reference=self.reference,
                            var_name="TCG761154TTA",
                            gene=self.aa2dna.get_gene("rpoB"),
                            mut="S450L")
        assert mutation.mutation_output_name == "S450L"

    def test_mutation_name_reverse_strand(self):
        gene = "gid"
        mutation_string = "I11N"
        is_protein_coding_var = True
        assert set(self.aa2dna.get_variant_names(
            gene, mutation_string, is_protein_coding_var)) == set(["GAT4408170ATT", "GAT4408170GTT"])
        mutation = Mutation(reference=self.reference,
                            var_name="GAT4408170ATT",
                            gene=self.aa2dna.get_gene("gid"),
                            mut="I11N")
        assert mutation.mutation_output_name == "I11N"

    def test_mutation_name_dna_space(self):
        gene = "pncA"
        mutation_string = "C18CCA"
        is_protein_coding_var = False
        assert set(self.aa2dna.get_variant_names(
            gene, mutation_string, is_protein_coding_var)) == set(["G2289224TGG"])
        mutation = Mutation(reference=self.reference,
                            var_name=self.aa2dna.get_variant_names(
                                gene, mutation_string, is_protein_coding_var)[0],
                            gene=self.aa2dna.get_gene(gene),
                            mut=mutation_string)
        assert mutation.mutation_output_name == "C18CCA"
