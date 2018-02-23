from __future__ import print_function
import sys
from collections import Counter
from mongoengine import connect
from mongoengine.connection import MongoEngineConnectionError
from pymongo.errors import ServerSelectionTimeoutError
import logging
from ga4ghmongo.schema import Variant
from mykatlas.utils import split_var_name
from mykatlas.utils import flatten
from mykatlas.utils import unique
from mykatlas.panelgeneration import AlleleGenerator
from ga4ghmongo.schema import VariantSet

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


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
    variant_to_samples = {}
    for variant in variants:
        variant_to_samples[variant] = variant.seen_in_samples()

    samples_counter = Counter(flatten(variant_to_samples.values()))
    samples_seen_more_than_once = [
        k for k,
        v in samples_counter.iteritems() if v > 1]
    contexts = []
    for sample in samples_seen_more_than_once:
        vars_together = []
        for variant, samples in variant_to_samples.items():
            if sample in samples:
                vars_together.append(variant)
        if vars_together not in contexts:
            contexts.append(vars_together)
            variants = [var for var in variants if var not in vars_together]
    for var in variants:
        contexts.append([var])
    return contexts + [[]]


def make_variant_probe(al, variant, kmer, DB=None, no_backgrounds=False):
    if no_backgrounds:
        context = []
    else:
        if DB is not None:
            try:
                context = get_context(variant.start, kmer)
            except (ServerSelectionTimeoutError, MongoEngineConnectionError):
                DB = None
                context = []
                logging.warning(
                    "Could not connect to database. Continuing without using genetic backgrounds")
        else:
            context = []
    if context:
        logging.info(
            "Found %i variants in context of %s" %
            (len(context), variant))
    variant_probe = None
    contexts_seen_together = seen_together(context)
    alts = []
    for context in contexts_seen_together:
        try:
            panel = al.create(variant, context)
        except ValueError as e:
            logging.warning("Failed to process variant:%s context:%s. %s" % (
                variant, ",".join([str(c) for c in context]), str(e)))
        else:
            ref = panel.ref
            panel.alts
            if variant_probe is not None:
                variant_probe.alts.extend(panel.alts)
            else:
                variant_probe = panel
    if variant_probe:
        variant_probe.alts = unique(variant_probe.alts)
    return variant_probe
