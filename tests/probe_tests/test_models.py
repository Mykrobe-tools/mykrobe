import sys
import json
import os

sys.path.append(".")
from mykrobe.probes.models import Mutation
from mykrobe.annotation.genes import GeneAminoAcidChangeToDNAVariants


class TestMutation():

    def setup(self):
        self.reference_filepath = "src/mykrobe/data/NC_000962.3.fasta"
        self.reference = os.path.basename(
            self.reference_filepath).split('.fa')[0]
        self.aa2dna = GeneAminoAcidChangeToDNAVariants(
            "src/mykrobe/data/NC_000962.3.fasta",
            "src/mykrobe/data/NC_000962.3.gb")

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
        assert mutation.mut == "S450L"

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
        assert mutation.mut == "I11N"
