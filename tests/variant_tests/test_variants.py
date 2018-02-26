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
        DB.drop_database('ga4ghmongo-test')

    def teardown(self):
        DB.drop_database('ga4ghmongo-test')


class TestVariantSets(BaseTest):

    def setup(self):
        DB.drop_database('ga4ghmongo-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")

    def test_create_new_variant(self):
        variant_set = VariantSet.create_and_save(
            name="this_vcf_file",
            reference_set=self.reference_set)
        vs = VariantSet.objects.get(name="this_vcf_file")
        assert variant_set == vs
        assert vs.reference_set.name == "ref_set"


class TestVariants(BaseTest):

    def setup(self):
        DB.drop_database('ga4ghmongo-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file",
            reference_set=self.reference_set)
        self.variant_sets = [self.variant_set]
        self.reference = Reference().create_and_save(
            name="ref",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])

    def teardown(self):
        DB.drop_database('ga4ghmongo-test')

    def test_create_SNP(self):
        v1 = Variant.create(variant_sets=self.variant_sets, start=0,
                            end=1, reference_bases="A",
                            alternate_bases=["T"],
                            reference=self.reference)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alternate_bases == ["T"]
        assert v1.length == 0

    def test_create_insertion(self):
        v1 = Variant.create(variant_sets=self.variant_sets,
                            start=0, end=1, reference_bases="T",
                            alternate_bases=["TA"],
                            reference=self.reference)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alternate_bases == ["TA"]
        assert v1.is_insertion
        assert v1.is_deletion is False
        assert v1.is_indel
        assert v1.length == 1

    def test_create_deletion(self):
        v1 = Variant.create(variant_sets=self.variant_sets,
                            start=0, end=1, reference_bases="AA",
                            alternate_bases=["A"],
                            reference=self.reference)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alternate_bases == ["A"]
        assert v1.reference_bases == "AA"
        assert v1.is_insertion is False
        assert v1.is_deletion
        assert v1.is_indel
        assert v1.length == 1

    def test_split_name(self):
        name = "A12T"
        r, pos, a = split_var_name(name)
        assert r == "A"
        assert pos == 12
        assert a == "T"

    def test_split_name_del(self):
        name = "AA12T"
        r, pos, a = split_var_name(name)
        assert r == "AA"
        assert pos == 12
        assert a == "T"

    def test_split_name_ins(self):
        name = "A12TT"
        r, pos, a = split_var_name(name)
        assert r == "A"
        assert pos == 12
        assert a == "TT"

    def test_split_name2(self):
        name = "A12T/A"
        r, pos, a = split_var_name(name)
        assert r == "A"
        assert pos == 12
        assert a == "T/A"

    def test_split_name3(self):
        name = "C-54T"
        r, pos, a = split_var_name(name)
        assert r == "C"
        assert pos == -54
        assert a == "T"
