import datetime
from mykrobe._vcf import VCF
from mykrobe.variants.schema.models import VariantSet
from mykrobe.variants.schema.models import VariantSetMetadata
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantCallSet
from mykrobe.variants.schema.models import VariantCall
from mongoengine import connect
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Reference

DB = connect('mykrobe-test')


class BaseTest():

    def setup(self):
        DB.drop_database('mykrobe-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.reference = Reference().create_and_save(
            name="NC_000962.3",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])

    def teardown(self):
        DB.drop_database('mykrobe-test')


class TestAddNewVariantSet(BaseTest):

    def test_add_new_vcf_variant_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        # We create a global variant set as well as one for the individual VCF
        assert VariantSet.objects().count() == 2
        vs = VariantSet.objects()[0]
        assert len(Variant.objects()[0].variant_sets) == 2
        assert vs.name == "test.vcf"


class TestAddNewVariantSetMetaData(BaseTest):

    def test_add_new_vcf_variant_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantSetMetadata.objects().count() >= 2
        assert VariantSetMetadata.objects(key="KMER").count() == 2


class TestAddNewCallSet(BaseTest):

    def test_add_new_call_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        # Only one callset but the callset should belong to multiple variant
        # sets
        assert VariantCallSet.objects().count() == 1
        assert VariantCallSet.objects()[
            0].created_at <= datetime.datetime.now()
        assert len(VariantCallSet.objects()[0].variant_sets) == 2


class TestVariantsAndCalls(BaseTest):

    def test_add_add_variants_and_calls(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantCall.objects().count() == 21
        assert Variant.objects().count() == 21


class TestAddSecondVCF(BaseTest):

    def setup(self):
        DB.drop_database('mykrobe-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.reference = Reference().create_and_save(
            name="NC_000962.3",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()

    def test_add_second_vcf_variant_set(self):
        # This VCF only has one Variant which is not in the first VCF
        vcf = VCF(
            f="tests/vcf_tests/test2.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantSet.objects().count() == 3
        assert VariantCallSet.objects().count() == 2
        assert VariantCall.objects().count() == 42
        assert Variant.objects().count() == 22
        assert len(Variant.objects()[0].variant_sets) == 3
        assert len(
            Variant.objects.get(
                names="UNION_BC_k31_var_147").variant_sets) == 3


class TestAddVCFwithIndels(BaseTest):

    def test_add_second_vcf_variant_set(self):
        # This VCF only has one Variant which is not in the first VCF
        vcf = VCF(
            f="tests/vcf_tests/test3.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantSet.objects().count() == 2
        assert VariantCallSet.objects().count() == 1
        assert VariantCall.objects().count() == 106
        assert Variant.objects().count() == 106
        assert Variant.snps().count() == 89
        assert Variant.indels().count() == 17
        assert Variant.insertions().count() == 8
        assert Variant.deletions().count() == 8
        assert Variant.ph_snps.count() == 1
