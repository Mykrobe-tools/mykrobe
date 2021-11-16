import logging

from mongoengine import DoesNotExist
from mongoengine import connect
from pymongo import MongoClient

from mykrobe._vcf import VCF
from mykrobe.constants import DB_PREFIX
from mykrobe.utils import get_first_chrom_name
from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet

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
    connect(DBNAME, host=args.db_uri)
    logger.debug("Using DB %s" % DBNAME)

    with open(args.reference_set) as fp:
        reference_set_name = get_first_chrom_name(fp)

    try:
        reference_set = ReferenceSet.objects.get(name=reference_set_name)
    except DoesNotExist:
        reference_set = ReferenceSet.create_and_save(name=reference_set_name)
        
    try:
        reference = Reference.create_and_save(
            name=reference_set_name, reference_sets=[reference_set], md5checksum="NA"
        )
    except:
        pass
    vcf = VCF(args.vcf, reference_set.id, method=args.method, force=args.force)
    vcf.add_to_database()
