from __future__ import print_function
import sys
from collections import Counter
from mongoengine import connect
from pymongo.errors import ServerSelectionTimeoutError
from mykrobe.variants.schema.models import Variant
from mykrobe.utils import split_var_name
from mykrobe.utils import flatten
from mykrobe.utils import unique
from mykrobe.probes import AlleleGenerator
from mykrobe.variants.schema.models import VariantSet

import logging
logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)


def get_context(pos, kmer):
    context = []
    for variant in Variant.objects(
            start__ne=pos,
            start__gt=pos - kmer,
            start__lt=pos + kmer):
        for split_variant in variant.split():
            context.append(split_variant)
    return context


def seen_together(variants):
    # Takes a list of variants.
    # Returns a list of variants that appear together (in the same variant set)
    variant_id_to_samples = {}
    for variant in variants:
        variant_id_to_samples[variant.var_hash] = variant.seen_in_samples()

    samples_counter = Counter(flatten(variant_id_to_samples.values()))
    samples_seen_more_than_once = [
        k for k,
        v in samples_counter.items() if v > 1]
    contexts = []

    for sample in samples_seen_more_than_once:
        vars_together = []

        for variant_id, samples in variant_id_to_samples.items():
            if sample in samples:
                vars_together.append(variant_id)

        if vars_together not in contexts:
            contexts.append(vars_together)
            variants = [
                var for var in variants if var.var_hash not in vars_together]

    for var in variants:
        contexts.append([var.var_hash])

    new_contexts = []
    for context in contexts:
        new_context = []
        for variant_id in context:
            new_context.append(Variant.objects.get(var_hash=variant_id))
        new_contexts.append(new_context)
    return new_contexts + [[]]


def make_variant_probe(al, variant, kmer, DB=None, no_backgrounds=False):
    if no_backgrounds:
        context = []
    else:
        if DB is not None:
            try:
                context = get_context(variant.start, kmer)
            except:
                DB = None
                context = []
                # logger.warning(
                # "Could not connect to database. Continuing without using backgrounds")
        else:
            context = []
    # if context:
        # logger.debug(
            # "Found %i variants in context of %s" %
            # (len(context), variant))
    variant_probe = None
    contexts_seen_together = seen_together(context)
    alts = []
    for context in contexts_seen_together:
        # logger.debug("Processing variant:%s context:%s" % (
            # variant, ",".join([str(c) for c in context])))
        try:
            panel = al.create(variant, context)
        except ValueError as e:
            pass
            # logger.warning("Failed to process variant:%s context:%s. %s" % (
            # variant, ",".join([str(c) for c in context]), str(e)))
        else:
            if variant_probe is not None:
                variant_probe.alts.extend(panel.alts)
                variant_probe.refs.extend(panel.refs)
            else:
                variant_probe = panel
    if variant_probe:
        variant_probe.alts = unique(variant_probe.alts)
        variant_probe.refs = unique(variant_probe.refs)
    return variant_probe
