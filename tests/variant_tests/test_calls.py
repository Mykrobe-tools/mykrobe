from ga4ghmongo.schema import Variant
from ga4ghmongo.schema import VariantSet
from ga4ghmongo.schema import VariantCall
from ga4ghmongo.schema import VariantCallSet
from ga4ghmongo.schema import Reference
from ga4ghmongo.schema import ReferenceSet

from ga4ghmongo.utils import split_var_name
from mongoengine import connect
DB = connect('ga4ghmongo-test')


class BaseTest():

    def setup(self):
        pass

    def teardown(self):
        DB.drop_database('ga4ghmongo-test')


class TestCallSet(BaseTest):

    def setup(self):
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file",
            reference_set=self.reference_set)
        self.variant_sets = [self.variant_set]

    def test_create_call_set(self):
        call_set = VariantCallSet.create_and_save(
            name="call_set",
            sample_id="C00123",
            variant_sets=self.variant_sets)
        cs = VariantCallSet.objects.get(name="call_set")
        assert call_set == cs
        assert cs.name == "call_set"
        assert cs.variant_sets[0].reference_set.name == "ref_set"


class TestCall(BaseTest):

    def setup(self):
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file2",
            reference_set=self.reference_set)
        self.variant_sets = [self.variant_set]
        self.reference = Reference().create_and_save(
            name="ref",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])
        self.call_set = VariantCallSet.create(
            sample_id="C00123",
            name="C00123",
            variant_sets=self.variant_sets)
        self.variant_snp = Variant.create(variant_sets=self.variant_sets,
                                          start=0, end=1, reference_bases="A",
                                          alternate_bases=["T"],
                                          reference=self.reference)

        self.variant_snp_mult_alts = Variant.create(
            variant_sets=self.variant_sets,
            start=0,
            end=1,
            reference_bases="T",
            alternate_bases=[
                "A",
                "C"],
            reference=self.reference)

    def test_create_SNP_het_call(self):
        c1 = VariantCall.create(variant=self.variant_snp,
                                call_set=self.call_set,
                                genotype=[0, 1],
                                genotype_likelihoods=[0.1, 0.9, 0.12])
        assert c1.call_set_name == "C00123"
        assert c1.genotype == [0, 1]
        self.variant_snp.save()
        c1.save()
        assert c1 in self.variant_snp.calls

        c2 = VariantCall.create(variant=self.variant_snp,
                                call_set=self.call_set,
                                genotype="1/1",
                                genotype_likelihoods=[0.01, 0.1, 0.9])
        assert c2.call_set_name == "C00123"
        assert c2.genotype == [1, 1]

        c2.save()
        assert c2 in self.variant_snp.calls

    def test_create_complex_call(self):
        c1 = VariantCall.create(
            variant=self.variant_snp_mult_alts,
            call_set=self.call_set,
            genotype="2/1",
            genotype_likelihoods=[
                0.01,
                0.1,
                0.9,
                0.1,
                0.2,
                0.6])
        self.variant_snp.save()
        self.variant_snp_mult_alts.save()
        c1.save()
        assert c1.call_set_name == "C00123"
        assert c1.genotype == [2, 1]
        assert c1 not in self.variant_snp.calls
        assert c1 in self.variant_snp_mult_alts.calls
