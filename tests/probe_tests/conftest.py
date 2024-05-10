import os

from mongoengine import connect
import pytest

from mykrobe.probes import AlleleGenerator
from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import VariantSet

DATA_DIR = os.path.join("tests", "ref_data")


@pytest.fixture(autouse=True)
def db_setup_teardown():
    DB = connect("mykrobe-test")
    DB.drop_database("mykrobe-test")
    yield
    DB.drop_database("mykrobe-test")


@pytest.fixture()
def variant_sets_and_reference():
    reference_set = ReferenceSet().create_and_save(name="ref_set")
    variant_set = VariantSet.create_and_save(
        name="this_vcf_file", reference_set=reference_set
    )
    variant_sets = [variant_set]
    reference = Reference().create_and_save(
        name="ref", md5checksum="sre", reference_sets=[reference_set]
    )
    return variant_sets, reference


@pytest.fixture()
def pg():
    pg = AlleleGenerator(reference_filepath=f"{DATA_DIR}/BX571856.1.fasta", kmer=31)
    return pg


@pytest.fixture()
def pg2():
    pg2 = AlleleGenerator(reference_filepath=f"{DATA_DIR}/NC_000962.3.fasta", kmer=31)
    return pg2
