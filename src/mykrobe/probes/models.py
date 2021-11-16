import os
from mongoengine import Document
from mongoengine import ReferenceField
from mongoengine import StringField
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from copy import copy
import itertools
from collections import Counter
import logging
import datetime
import math

from mykrobe import K
from mykrobe.utils import make_hash
from mykrobe.utils import split_var_name
from mykrobe.utils import seq_to_kmers

from mykrobe.variants.schema.models.base import CreateAndSaveMixin
from mykrobe.variants.schema.models import Variant


def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]


class VariantPanel(Document, CreateAndSaveMixin):

    meta = {'indexes': [
        {
            'fields': ['variant']
        },
        {
            'fields': ['var_hash']
        }

    ]
    }

    variant = ReferenceField('Variant')
    var_hash = StringField()

    @classmethod
    def create(cls, variant):
        return cls(
            variant=variant,
            var_hash=variant.var_hash
        )

    def __repr__(self):
        return "VariantPanel %s" % self._name


class AlleleGenerator(object):

    """docstring for PanelGenerator"""

    def __init__(self, reference_filepath, kmer=K):
        self.reference_filepath = reference_filepath
        self.kmer = kmer
        self.ref = []
        self.ref_length = 0
        self._read_reference()

    def _read_reference(self):
        for record in SeqIO.parse(self.reference_filepath, 'fasta'):
            self.reference_sequence=record.seq.upper()
            self.ref += list(record.seq.upper())
        ### Pad with Ns for SNPs at the end of the reference
        self.ref.extend(["N"]*self.kmer)
        self.ref_length = len(self.ref)

    def create(self, v, context=[]):
        # Position should be 1 based
        self._check_valid_variant(v)
        context = self._remove_overlapping_contexts(v, context)
        context = self._remove_contexts_not_within_k(v, context)
        wild_type_reference = self._get_wildtype_reference(v)
        null_variant = Variant.create(
            v.start, v.reference_bases, alternate_bases=[v.reference_bases])
        references = self._generate_alternates_on_all_backgrounds(
            null_variant, context)
        alternates = self._generate_alternates_on_all_backgrounds(v, context)
        ## The alternates shouldn't contain kmers in the reference
        if v.is_indel:
            alternates=self.trim_uninformative_kmers(alternates, references)
        return Panel(v, references, v.start, alternates)

    def trim_uninformative_kmers(self, alternates, references=[]):
        new_alternates=[]
        for ref, alt in zip(references, alternates):
            alt="".join(alt)
            ref="".join(ref)
            informative_kmers=[]
            for i,k in enumerate(seq_to_kmers(alt, self.kmer)):
                if not k in ref:
                    informative_kmers.append(i)
            if informative_kmers:
                trim=(min(informative_kmers), max(informative_kmers))
                alt=alt[trim[0]:trim[1]+self.kmer]
            assert alt

            informative_kmers=[]
            for i,k in enumerate(seq_to_kmers(alt, self.kmer)):
                if not k in self.reference_sequence:
                    informative_kmers.append(i)
            if informative_kmers:
                trim=(min(informative_kmers), max(informative_kmers))
                alt=alt[trim[0]:trim[1]+self.kmer]
            assert alt

            new_alternates.append(alt)
        return new_alternates




    def _check_valid_variant(self, v):
        index = v.start - 1
        if not len(v.alternate_bases) == 1:
            raise NotImplementedError(
                "Probes can only be built for homozygous variants at this time")
        if "".join(
                self.ref[index:(index + len(v.reference_bases))]) != v.reference_bases:
            raise ValueError(
                "Cannot create alleles as ref at pos %i is not %s (it's %s) are you sure you're using one-based co-ordinates?" %
                (v.start,
                 v.reference_bases,
                 "".join(
                     self.ref[
                         index: (
                             index +
                             len(
                                 v.reference_bases))])))
        if v.start <= 0:
            raise ValueError("Position should be 1 based")

    def _remove_overlapping_contexts(self, v, context):
        """If there's a var within the range of an DEL remove from context"""
        context = [c for c in context if not c.overlapping(v)]
        return context

    def _remove_contexts_not_within_k(self, v, context):
        new_context = []
        for c in context:
            if c.is_insertion:
                effective_pos = c.start - c.length
            elif c.is_deletion:
                effective_pos = c.start + c.length
            else:
                effective_pos = c.start
            if abs(v.start - effective_pos) < self.kmer:
                new_context.append(c)
        return new_context

    def _get_wildtype_reference(self, v):
        i, start_index, end_index = self._get_start_end(v)
        return self._get_reference_segment(start_index, end_index)

    def _get_reference_segment(self, start_index, end_index):
        return copy(self.ref[start_index:end_index])

    def _get_alternate_reference_segment(self, v, context):
        ref_segment_length_delta = self._calculate_length_delta_from_indels(v,            context)
        i, start_index, end_index = self._get_start_end(
            v, delta=ref_segment_length_delta)
        return self._get_reference_segment(start_index, end_index)

    def _generate_alternates_on_all_backgrounds(self, v, context):
        # First, create all the context combinations
        context_combinations = self._get_all_context_combinations(context)
        # For each context, create the background and alternate
        alternates = []
        for context_combo in context_combinations:
            ref_segment_length_delta = self._calculate_length_delta_from_indels(v, context_combo)
            i, start_index, end_index = self._get_start_end(
                v, delta=ref_segment_length_delta)
            alternate_reference_segment = self._get_reference_segment(
                start_index,
                end_index)
            try:
                background = self._generate_background_using_context(
                    i,
                    v,
                    alternate_reference_segment,
                    context_combo)
            except (ValueError, AssertionError) as e:
                m = "Could not process context combo %s. " % (
                    ",".join([c.var_name for c in context_combo] + [v.var_name]))
                raise ValueError("\n".join([m, str(e)]))
            alternate = copy(background)
            i -= self._calculate_length_delta_from_variant_list(
                [c for c in context_combo if c.start <= v.start and c.is_indel])
            if not "".join(
                    alternate[i:(i + len(v.reference_bases))]) == v.reference_bases:
                raise ValueError("Could not process context combo %s. %s != %s " %
                                 (",".join([c.var_name for c in context_combo] +
                                           [v.var_name]), "".join(alternate[i:(i +
                                                                               len(v.reference_bases))]), v.reference_bases))
            else:
                for alt in v.alternate_bases:
                    alternate[i: i + len(v.reference_bases)] = alt
                    alternates.append(alternate)
        return alternates

    def _get_all_context_combinations(self, context):
        context_list = [[]]
        if context:
            # Create contexts if multiple variants at given posision
            contexts = self._create_multiple_contexts(context)
            for context in contexts:
                context_combinations = self._get_combinations_of_backgrounds(
                    context)
                context_list.extend(context_combinations)
        return context_list

    def _generate_background_using_context(
            self,
            i,
            v,
            alternate_reference_segment,
            context):
        backgrounds = [alternate_reference_segment]
        new_background = copy(alternate_reference_segment)

        variants_added = []
        for e, variant in enumerate(context):
            for alt in variant.alternate_bases:
                j = i + variant.start - v.start
                j -= self._calculate_length_delta_from_variant_list(
                    [c for c in variants_added if c.start <= variant.start and c.is_indel])
                if j <= len(new_background) and j >= 0:
                    if j + len(variant.reference_bases) > len(new_background):
                        hang = j + \
                            len(variant.reference_bases) - len(new_background)
                        assert "".join(
                            new_background[
                                j: len(new_background)]) == variant.reference_bases[
                            :len(
                                variant.reference_bases) -
                            hang]
                        new_background[
                            j: j +
                            len(
                                variant.reference_bases)] = alt[
                            :len(
                                variant.reference_bases) -
                            hang]
                    else:
                        if not "".join(new_background[
                                j: j + len(variant.reference_bases)]) == variant.reference_bases:
                            raise ValueError(
                                "Could not process variant %s. %s != %s " %
                                (variant.var_name,
                                 "".join(
                                     new_background[
                                         j: j +
                                         len(
                                             variant.reference_bases)]),
                                    variant.reference_bases))
                        else:
                            new_background[
                                j: j + len(variant.reference_bases)] = alt
                    variants_added.append(variant)
                else:
                    del context[e]
        return new_background

    def _create_multiple_contexts(self, context):
        new_contexts = self._recursive_context_creator([context])
        return new_contexts

    def _recursive_context_creator(self, contexts):
        # This is only run when there are multiple variants at the same
        # position
        compatiblity_of_contexts = [
            self._all_variants_are_combatible(context) for context in contexts]
        if self._are_contexts_all_valid(compatiblity_of_contexts):
            return contexts
        else:
            # If there are variants that are incompatible (e.g. SNPs at the same position or overlapping INDELs)
            # we need to distangle them into multiple contexts which are all compatible with each other
            # Get the first incompatible context
            i = compatiblity_of_contexts.index(False)
            incompatible_context = contexts.pop(i)
            # Split it into two smaller contexts by removing incompatible
            # variants
            split_context = self._split_context(incompatible_context)
            contexts.extend(split_context)
            return self._recursive_context_creator(contexts)

    def _split_context(self, context):
        assert not self._all_variants_are_combatible(context)
        v1, v2 = self._get_first_two_incompatible_variants(context)
        nov1 = []
        nov2 = []
        for x in context:
            if x is not v1:
                nov1.append(x)
            if x is not v2:
                nov2.append(x)
        return [nov1, nov2]

    def _get_first_two_incompatible_variants(self, context):
        for pair in itertools.combinations(context, 2):
            if pair[0].overlapping(pair[1]):
                return pair

    def _are_contexts_all_valid(self, compatiblity_of_contexts):
        return all(compatiblity_of_contexts)

    def _all_variants_are_combatible(self, variants):
        pairs = list(itertools.combinations(variants, 2))
        return all([not v1.overlapping(v2) for v1, v2 in pairs])

    def _get_combinations_of_backgrounds(self, context):
        combination_context = []
        for L in range(1, len(context) + 1):
            for subset in itertools.combinations(context, L):
                combination_context.append(list(subset))
        return combination_context

    def _get_start_end(self, v, delta=0):
        # Is large var
        shift = 0
        kmer = self.kmer
        if len(v.reference_bases) > 2 * kmer:
            kmer = int(math.ceil(float(len(v.reference_bases)) / 2)) + 5
        elif (v.length > 2 * kmer):
            kmer = int(math.ceil(float(v.length) / 2)) + 5
        if len(v.reference_bases) > kmer:
            shift = int(
                (kmer - 1) - math.floor(float((2 * kmer + 1) - len(v.reference_bases)) / 2))
        pos = v.start
        start_delta = int(math.floor(float(delta) / 2))
        end_delta = int(math.ceil(float(delta) / 2))
        start_index = pos - kmer - start_delta
        end_index = pos + kmer + end_delta-1
        min_probe_length=(2 * kmer) - 1
        i = kmer - 1 + start_delta
        ### Is the variant at the start of the sequence? This is a special case.
        if start_index < 0:
            diff = abs(start_index)
            start_index = 0
            end_index += diff
            i -= diff
        start_index += shift
        end_index += shift
        i -= shift
        # print(start_index, end_index)
        if (end_index - start_index) >= min_probe_length:
            return (i, start_index, end_index)
        else:
            return self._get_start_end(v, delta=0)

    def _calculate_length_delta_from_indels(self, v, context):
        """Calculates the change in required bases for given variant.
        For deletions we need extra bases to get same length flanks"""
        variants = context + [v]
        return self._calculate_length_delta_from_variant_list(variants)

    def _calculate_length_delta_from_variant_list(self, variants):
        deletions_length = [v.length for v in variants if v.is_deletion]
        insertions_length = [v.length for v in variants if v.is_insertion]
        delta = sum(deletions_length) - sum(insertions_length)
        if delta < -self.kmer:
            return abs(delta)
        else:
            return delta


class Panel(object):

    def __init__(self, variant, refs, start, alts):
        self.variant = variant
        self.refs = unique(["".join(ref) for ref in refs])
        self.start = start
        self.alts = unique(["".join(alt) for alt in alts])
        self.alts=list(set(self.alts)-set(self.refs))


class Mutation(object):

    def __init__(
            self,
            var_name,
            reference,
            gene=None,
            mut=None,
            protein_coding_var=False):
        self.var_name = var_name
        self.gene = gene
        if mut:
            tmp, self.start, tmp = split_var_name(mut)
        self.ref, tmp, self.alt = split_var_name(var_name)
        self.standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.reference = reference
        self.input_mutation_name = mut
        self.protein_coding_var = protein_coding_var

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    @property
    def mutation_output_name(self):
        if self.input_mutation_name:
            return self.input_mutation_name
        else:
            if self.gene is not None:
                if self.gene.forward:
                    ref = self.ref
                    alt = self.alt
                else:
                    ref = str(Seq(self.ref).reverse_complement())
                    alt = str(Seq(self.alt).reverse_complement())
                if self.protein_coding_var:
                    r = self.standard_table.forward_table.get(ref, ref)
                    a = self.standard_table.forward_table.get(alt, alt)
                else:
                    r = self.ref
                    a = self.alt
                return "".join([r, str(self.start), a])
            else:
                return self.var_name

    @property
    def variant(self):
        ref, start, alt = split_var_name(self.var_name)
        return Variant.create(variant_sets=None, start=int(start),
                              end=0, reference_bases=ref,
                              alternate_bases=[alt],
                              reference=self.reference)
