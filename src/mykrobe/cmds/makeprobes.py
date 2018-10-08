from __future__ import print_function
import csv
import os
import sys
import logging
from mongoengine import connect

from mongoengine import DoesNotExist
from mongoengine import NotUniqueError

from pymongo.errors import ServerSelectionTimeoutError
from Bio.Seq import Seq


from atlas.variants.schema.models import Variant
from atlas.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Reference

from mykrobe.utils import split_var_name
from mykrobe.annotation.genes import GeneAminoAcidChangeToDNAVariants
from mykrobe._vcf import VCF

from mykrobe.probes.models import Mutation
from mykrobe.probes import AlleleGenerator
from mykrobe.probes import make_variant_probe
from mykrobe.constants import DB_PREFIX


logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


def run(parser, args):
    DB = connect('%s-%s' % (DB_PREFIX, args.db_name))
    if DB is not None:
        try:
            Variant.objects()
            logging.info(
                "Connected to %s-%s" % (DB_PREFIX, args.db_name))
        except (ServerSelectionTimeoutError):
            DB = None
            logging.warning(
                "Could not connect to database. Continuing without using genetic backgrounds")
    mutations = []
    reference = os.path.basename(args.reference_filepath).split('.fa')[0]
    if args.vcf:
        run_make_probes_from_vcf_file(args)
    elif args.genbank:
        aa2dna = GeneAminoAcidChangeToDNAVariants(
            args.reference_filepath,
            args.genbank)
        if args.text_file:
            with open(args.text_file, 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                for row in reader:
                    gene, mutation_string, alphabet = row
                    if alphabet == "DNA":
                        protein_coding_var = False
                    else:
                        protein_coding_var = True
                    for var_name in aa2dna.get_variant_names(
                            gene, mutation_string, protein_coding_var):
                        mutation = Mutation(reference=reference,
                                            var_name=var_name,
                                            gene=aa2dna.get_gene(gene),
                                            mut=mutation_string,
                                            protein_coding_var=protein_coding_var)
                        mutations.append(mutation)
        else:
            for variant in args.variants:

                gene, mutation = variant.split("_")
                for var_name in aa2dna.get_variant_names(gene, mutation):
                    mutations.append(
                        Mutation(reference=reference,
                                 var_name=var_name,
                                 gene=gene,
                                 mut=mutation))
    else:
        if args.text_file:
            with open(args.text_file, 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                for row in reader:
                    gene_name, pos, ref, alt, alphabet = row
                    if gene_name == "ref":
                        mutations.append(
                            Mutation(
                                reference=reference,
                                var_name="".join([ref, pos, alt])))
                    else:
                        mutations.append(
                            Mutation(
                                reference=reference,
                                var_name=row[0]))
        else:
            mutations.extend(Mutation(reference=reference, var_name=v)
                             for v in args.variants)
    al = AlleleGenerator(
        reference_filepath=args.reference_filepath,
        kmer=args.kmer)
    for enum, mut in enumerate(mutations):
        if enum % 100 == 0:
            logger.info(
                "%i of %i - %f%%" % (enum, len(mutations), round(100*enum/len(mutations), 2)))
        variant_panel = make_variant_probe(
            al, mut.variant, args.kmer, DB=DB, no_backgrounds=args.no_backgrounds)
        if variant_panel is not None:
            for i, ref in enumerate(variant_panel.refs):
                try:
                    gene_name = mut.gene.name
                except AttributeError:
                    gene_name = "NA"

                sys.stdout.write(
                    ">ref-%s?var_name=%s&num_alts=%i&ref=%s&enum=%i&gene=%s&mut=%s\n" %
                    (mut.mutation_output_name, mut.variant.var_name, len(
                        variant_panel.alts), mut.reference, i, gene_name, mut.mutation_output_name))
                sys.stdout.write("%s\n" % ref)

            for i, a in enumerate(variant_panel.alts):
                sys.stdout.write(">alt-%s?var_name=%s&enum=%i&gene=%s&mut=%s\n" %
                                 (mut.mutation_output_name, mut.variant.var_name, i, gene_name, mut.mutation_output_name))

                sys.stdout.write("%s\n" % a)
        else:
            logging.warning(
                "All variants failed for %s_%s - %s" %
                (mut.gene, mut.mutation_output_name, mut.variant))


def run_make_probes_from_vcf_file(args):
    # Make VariantSet from vcf
    reference = os.path.basename(args.reference_filepath).split(".fa")[0]
    try:
        reference_set = ReferenceSet.objects.get(name=reference)
    except DoesNotExist:
        reference_set = ReferenceSet.create_and_save(name=reference)
        # Hack
    try:
        reference = Reference.create_and_save(
            name=reference,
            reference_sets=[reference_set],
            md5checksum=reference)
    except NotUniqueError:
        pass
    vcf = VCF(
        args.vcf,
        reference_set.id,
        method="tmp",
        force=True,
        append_to_global_variant_set=False)
    vcf.add_to_database()
