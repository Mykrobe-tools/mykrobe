from __future__ import print_function

import itertools
import logging
from collections import defaultdict

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

from mykrobe import STOP
from mykrobe.utils import split_var_name

logger = logging.getLogger(__name__)
logger.setLevel("INFO")


def flatten(l):
    return [item for sublist in l for item in sublist]


class Region(object):
    def __init__(self, reference, start, end, forward=True):
        self.reference = reference
        self.start = start
        self.end = end
        self.forward = forward

    @property
    def strand(self):
        if self.forward:
            return "forward"
        else:
            return "reverse"

    @property
    def seq(self):
        if self.forward:
            return self.reference[self.start - 1 : self.end]
        else:
            return self.reference[self.start - 1 : self.end].reverse_complement()

    def get_reference_position(self, pos):
        if pos < 0 and self.forward:
            return self.start + pos
        elif pos < 0 and not self.forward:
            # Upstream of a gene on the reverse stand
            return self.end - pos
        elif pos > 0 and self.forward:
            return self.start + pos - 1
        elif pos > 0 and not self.forward:
            return self.end - pos + 1
        else:
            raise ValueError("Positions are 1-based")


class Gene(Region):
    def __init__(self, name, reference, start, end, forward=True):
        super(self.__class__, self).__init__(reference, start, end, forward)
        self.name = name
        self.translation_table = 11
        self.backward_codon_table = make_backward_codon_table()

    @property
    def prot(self):
        return self.seq.translate(table=self.translation_table).rstrip("*")

    def get_context(self, pos, N):
        return self.seq[(3 * (pos - 1)) - N : pos * 3 + N]

    def get_codon(self, pos):
        if pos > len(self.prot):
            raise ValueError(
                "There are only %s aminoacids in this gene" % len(self.prot)
            )
        else:
            return self.seq[(3 * (pos - 1)) : pos * 3]

    def get_reference_codon(self, pos):
        if self.forward:
            return self.get_codon(pos)
        else:
            return self.get_codon(pos).reverse_complement()

    def get_reference_codons(self, pos):
        standard_table = CodonTable.unambiguous_dna_by_id[11]
        # Get the reference codon in reading frame
        ref_codon = self.get_codon(pos)
        # Get the backward codons
        ref_aa = standard_table.forward_table[str(ref_codon)]
        condons_in_reading_frame = self.backward_codon_table[ref_aa]
        if self.forward:
            return condons_in_reading_frame
        else:
            return [str(Seq(s).reverse_complement()) for s in condons_in_reading_frame]

    def __str__(self):
        return "Gene:%s" % self.name

    def __repr__(self):
        return "Gene:%s" % self.name


def make_backward_codon_table():
    table = defaultdict(list)
    standard_table = CodonTable.unambiguous_dna_by_id[11]
    codons = generate_all_possible_codons()
    for codon in codons:
        if codon not in standard_table.stop_codons:
            try:
                table[standard_table.forward_table[codon]].append(codon)
            except:
                table[standard_table.forward_table[codon]] = [codon]
        else:
            table["*"].append(codon)
    return table


def generate_all_possible_codons():
    return [
        "".join(subset) for subset in itertools.product(["A", "T", "C", "G"], repeat=3)
    ]


class GeneAminoAcidChangeToDNAVariants:
    def __init__(self, reference, genbank):
        self.reference = self._parse_reference(reference)
        self.genbank = self._parse_genbank(genbank)
        self.backward_codon_table = make_backward_codon_table()

    def _parse_reference(self, reference):
        with open(reference, "r") as infile:
            return list(SeqIO.parse(infile, "fasta"))[0].seq.upper()

    def _parse_genbank(self, genbank):
        d = {}
        with open(genbank, "r") as infile:
            for feat in SeqIO.read(infile, "genbank").features:
                if feat.type == "gene":
                    try:
                        name = feat.qualifiers.get(
                            "gene", feat.qualifiers.get("db_xref")
                        )[0]
                    except TypeError:
                        pass
                    else:
                        # SeqIO converts to 0-based
                        strand = int(feat.location.strand)
                        if strand == 1:
                            forward = True
                        elif strand == -1:
                            forward = False
                        else:
                            raise ValueError("Strand must be 1 or -1")
                        d[name] = Gene(
                            name,
                            self.reference,
                            start=feat.location.start + 1,
                            end=feat.location.end,
                            forward=forward,
                        )
        return d

    def get_alts(self, amino_acid):
        if amino_acid == "X":
            # exclude stop codons
            non_stop_codons = [
                v for k, v in self.backward_codon_table.items() if k != STOP
            ]
            return flatten(non_stop_codons)
        else:
            return self.backward_codon_table[amino_acid]

    def get_reference_alts(self, gene, amino_acid):
        if gene.forward:
            return self.get_alts(amino_acid)
        else:
            return [str(Seq(s).reverse_complement()) for s in self.get_alts(amino_acid)]

    def get_location(self, gene, pos):
        if gene.forward:
            dna_pos = (3 * (pos - 1)) + 1
        else:
            dna_pos = 3 * (pos)
        return gene.get_reference_position(dna_pos)

    def get_variant_names(self, gene, mutation, protein_coding_var=True):
        ref, start, alt = split_var_name(mutation)
        gene = self.get_gene(gene)
        if start < 0 or not protein_coding_var:
            return self._process_DNA_mutation(gene, ref, start, alt)
        elif start > 0:
            return self._process_coding_mutation(gene, ref, start, alt)
        else:
            raise ValueError(
                "Variants are defined in 1-based coordinates. You can't have pos 0. "
            )

    def _process_DNA_mutation(self, gene, ref, start, alt):
        names = []
        pos = gene.get_reference_position(start)
        if not gene.forward:
            pos -= len(ref) - 1
            ref = str(Seq(ref).reverse_complement())
            if alt != "X":
                alt = str(Seq(alt).reverse_complement())
        if alt == "X":
            for a in ["A", "T", "C", "G"]:
                if a != ref:
                    names.append("".join([ref, str(pos), a]))
        else:
            names.append("".join([ref, str(pos), alt]))
        return names

    def _process_coding_mutation(self, gene, ref, start, alt):
        logger.debug(
            "Processing gene:{} ref:{} start:{} alt:{}".format(gene, ref, start, alt)
        )
        if not gene.prot or start > len(gene.prot):
            raise ValueError(
                "Error translating %s_%s " % (gene, "".join([ref, str(start), alt]))
            )
        if not gene.prot[start - 1] == ref:
            raise ValueError(
                "Error processing %s_%s. The reference at pos %i is not %s, it's %s. "
                % (
                    gene,
                    "".join([ref, str(start), alt]),
                    start,
                    ref,
                    gene.prot[start - 1],
                )
            )
        ref_codons = gene.get_reference_codons(start)
        alt_codons = self.get_reference_alts(gene, alt)
        logger.debug(
            "Reference codons (forward strand equivalent) {}".format(
                "".join(ref_codons)
            )
        )
        logger.debug(
            "Alternate codons (forward strand equivalent) {}".format(
                "".join(alt_codons)
            )
        )
        for ref_codon in ref_codons:
            if ref_codon in alt_codons:
                alt_codons.remove(ref_codon)
        location = self.get_location(gene, start)
        alternative = "/".join(alt_codons)
        ref_codon = gene.get_reference_codon(start)
        names = [
            "".join(["".join(ref_codon), str(location), "".join(alt_codon)])
            for alt_codon in alt_codons
        ]
        return names

    def get_gene(self, gene):
        return self.genbank[gene]
