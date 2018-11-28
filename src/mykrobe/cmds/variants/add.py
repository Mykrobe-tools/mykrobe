import logging
import os
import csv
from mongoengine import connect
from mongoengine import NotUniqueError
from mongoengine import OperationError
from mongoengine import DoesNotExist
from pymongo import MongoClient

from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Reference
from mykrobe._vcf import VCF
from mykrobe.constants import DB_PREFIX

"""Adds variants to the database"""

logger = logging.getLogger(__name__)
client = MongoClient()


def is_record_valid(record):
    valid = True
    for sample in record.samples:
        if sample["GT"] is None:
            valid = False
        else:
            if sum([int(i) for i in sample["GT"].split("/")]) < 2:
                valid = False
        try:
            if sample["GT_CONF"] <= 1:
                valid = False
        except AttributeError:
            pass
    return valid


def get_genotype_likelihood(sample):
    try:
        genotype_likelihood = sample["GT_CONF"]
    except AttributeError:
        genotype_likelihood = sample["GQ"]
    return genotype_likelihood


def run(parser, args):
    args = parser.parse_args()
    # args = check_args(args)
    if args.quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    DBNAME = "%s-%s" % (DB_PREFIX, args.db_name)
    db = client[DBNAME]
    connect(DBNAME)
    logger.debug("Using DB %s" % DBNAME)
    reference_set_name = ".".join(os.path.basename(args.reference_set).split(".")[:-1])
    try:
        reference_set = ReferenceSet.objects.get(name=reference_set_name)
    except DoesNotExist:
        reference_set = ReferenceSet.create_and_save(name=reference_set_name)
        # Hack
    try:
        reference = Reference.create_and_save(
            name=reference_set_name, reference_sets=[reference_set], md5checksum="NA"
        )
    except:
        pass
    vcf = VCF(args.vcf, reference_set.id, method=args.method, force=args.force)
    vcf.add_to_database()
