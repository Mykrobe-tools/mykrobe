from mongoengine import connect
import pytest

from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import VariantSet


@pytest.fixture(autouse=True)
def db_setup_teardown():
    DB = connect("mykrobe-test")
    DB.drop_database("mykrobe-test")
    yield
    DB.drop_database("mykrobe-test")


@pytest.fixture()
def reference_set():
    reference_set = ReferenceSet().create_and_save(name="ref_set")
    return reference_set


@pytest.fixture()
def reference(reference_set):
    return Reference().create_and_save(
        name="ref", md5checksum="sre", reference_sets=[reference_set]
    )


@pytest.fixture()
def variant_sets(reference_set, reference):
    variant_set = VariantSet.create_and_save(
        name="this_vcf_file", reference_set=reference_set
    )
    return [variant_set]
