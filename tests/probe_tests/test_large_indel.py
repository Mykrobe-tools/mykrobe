from mykrobe.probes import AlleleGenerator
from mykrobe.variants.schema import Variant
from mykrobe.variants.schema import VariantSet
from mykrobe.variants.schema import Reference
from mykrobe.variants.schema import ReferenceSet
from nose.tools import assert_raises
from mongoengine import connect
DB = connect('atlas-test')


class TestLargeINDELAlleleGenerator():

    def setup(self):
        DB.drop_database('atlas-test')
        self.pg = AlleleGenerator(
            reference_filepath="atlasvar/data/NC_000962.2.fasta")
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

    def test_large_variant(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACT",
            start=1355983,
            alternate_bases=["ACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCG"])
        panel = self.pg.create(v)
        assert panel.ref == "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"
        assert panel.alts == [
            "TCGTCACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCGGGACC"]

    def test_large_variant2(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="AACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCAC",
            start=1355983,
            alternate_bases=["ACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGC"])
        panel = self.pg.create(v)
        assert panel.ref == "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"
        assert panel.alts == [
            "TCGTCACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCTGGACC"]

    def test_large_variant3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC",
            start=1355978,
            alternate_bases=["TCGTCAACGCCCGGTATCTGAGGATCGGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"])
        panel = self.pg.create(v)
        assert panel.ref == "CAAACCTCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACCGCCATA"
        assert panel.alts == [
            "CAAACCTCGTCAACGCCCGGTATCTGAGGATCGGTGTTCTCACCCAATACAAGTCGCATTCACTGGACCGCCATA"]

    def test_very_large_variant3(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACCGCCAT",
            start=1355978,
            alternate_bases=["TCGTCAACGCCCGGTATCTGAGGATCGGTGTTCACCCAATACAAGTCGCATTCACTGGACCGCCAT"])
        panel = self.pg.create(v)
        assert panel.ref == "AAACCTCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACCGCCATATCTCG"
        assert panel.alts == [
            "CAAACCTCGTCAACGCCCGGTATCTGAGGATCGGTGTTCACCCAATACAAGTCGCATTCACTGGACCGCCATATCTCGC"]

    def test_large_insertion(self):
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference,
            reference_bases="C",
            start=2352065,
            alternate_bases=["CCTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTG"])
        panel = self.pg.create(v)
        assert panel.ref == "ACCGGTCCAGCTCGGCCAGCTCAGTCACGTCGCCGCCGCCTCGCCAGTTGACCGCGCCCGCTCGCGGCTAGCGGGCCTA"
        assert panel.alts == [
            "AAGTCGTCGAGCGAGAACGGTAGTTCCGCGGTGAACCGGTCCAGCTCGGCCAGCTCAGTCACGTCGCCGCCGCCTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTGCTCGCCAGTTGACCGCGCCCGCTCGCGGCTAGCGGGCCTACGTGACGTCGTCATGAGATCCGATGACCGATGGC"]

    def test_large_var1(self):
        v = Variant.create(variant_sets=self.variant_sets, reference=self.reference, reference_bases="CGCGGGAGTAGAACGATCGCCAAGTGGTCGGTCTTGGCTGCCCACTTCATCCCCGGCGCCACCGGCAGGTCTCGCGGTCATCTCGACCAACGGAGGGCCGTCGGTGGTTCGTATCCGGCCAAGAACGGCGAGAACGGTTTGTGCCTCTATGCCAGGGTGAATGTCTCATCTCCCAGGCGGACGGTGATATCCAGTTCTCCGCCAAGAGCGGACACGTATTTGCGCAGTGTGTTGACCTGTGCGGAGCCGATGTCGCCGTTCTCGATGCTGGATACCCGGCTCTGCCGGATGTGCGCCAGCGCAGCCACCTGGACCTGGGTGAGTGACTGAGCCGCGCGCAGCTCCCGGAGCCGGAATGCCCGCACTTCATCGCGCATTCGTGCCTTGTGCCGGTCCACCGCCTCCCGGTTAACGGGACGTACGGCGTCCATGTCCCGTAGTGTCATCGCCATCGTGCCACTTACCCTTTCTTGCGCTTGCGCCTCTTTGGCTTCGTGTCCTCGAACTGTGCGAGATGTTCGGCAAACATCTCATCGGCCGCTTTGATCTTCTCGTCGTACCACTGGGTCCACCGCCCGGCCTTGTTACCGGCGGCCAGCATGATCGCCTGCCGCGCCGGGTCGAAGGCGAACAGAATGCGGACCTCGGACCGCCCTTGTGATCCTGGACGCAGCTCCTTCATGTTCTTGTGGCGCGACCCACGCACCGTGTCCACCAGAGGACAGCCAAGTGCGGGGCCCTCTTCCTCGAGAACCTCGATAGCTGCGAACACCAATTCGTAGGTCTCTCGGTCCAAGCCGTTGAGCCAGGCGGAGATGCGCTCCACATCCGCCGTCCACCCCACAGAGTCGCAGAGTAGCGCGATACGCGATATCACACAAGGGTGATATTCCTCCGGGTAAGAGCAGCGGGCGACGGGGCTACCGTCGAGGAAATGCCGGCAGGCGAGGACGGACTCTGCGCACCCGGGCCGTTGAAACAGTAGCCTGTGCCAGGCCGAGAATTCATCCCCACGTATGAGGCAGTACAGTGCGCCGCCGTGCGCGTTCTCCCATGGAACGTTCACGGGCTCCCGTGGATGACAGGCGTTTCATGAACGCCAGCGCCGCCGCAACCCGACCGAAAGCGGTTGACCCCAAGGAGAGCTGGAAGTCGAGGCCACCACCTTCGCCGCGGAGTTGCTCATGCCCGAGAGCGAGACTCGTCCCGAAATACGCCGGCTCGATTTCGGCAAGTTGCTCGAACTGAAGCGGGAATGGGCGTCGACCCGCTCGACCAGCCCCAGCCGGGTGACCAGCCCCAGCCGGGTGACCAGCCGATGCACCGCGGCGATCCCACCGAAGCCGGTGGCATCGATGTTGGCGCCGACCTCGTAGCGCACCGCGCCCGAACCCAGCATCGGCCTGGGCTGCGCCGCCCAGCGTCCAGCCCGCGCGTGCCGCGCCGCCACCCTGCGCCCTCGGCGTGTGATGTTTCGCCGACTCTGTTCATGGGTTATCTTCTTCACCACAAAGGCCTTTCCTGCTGGGCTGTGTTGAGGTCGCAAACCCAGCCAGGGTAAGGCCTTTGGCCTCTCCTACCCGGCCGACACGCTTACTGAAGGCCTAGTCTAGGCAGGCCATTCAATCTGCGGAATCGAAAAATTCGGTTCCAGCCTGCTCGTTTCCTTTCCGACAGCGATCTGACGTTGCGTAACGTCATTTGTACGGACTCTTTTAGCGGCATTGATTTCAGATGCCAACGCCGTCTGTGCTGTAGCGCCGATTGGCCGAAACTGTAAATTTGTATGATTATTTAAATCTTTGACGAACACGCGCCACAAACGTACTATCTCTTTGGCAAAGTCCACCGGCATCTCATTCAACGGTTTTGTTTGCGCGTGGTCGTCATATGTTGGTAACTGTGTAACCGGCCGCCTATCTTGCGCGTGCATCATATGACTATGAATCGGCCTTCTCCAGTGAAATTGATACAAGATCGATCCGATAAGCGGTACCTTGTACACAGTGCAATTGTAGTAATTCGCGTTTTGTCCTACGCTTGTATTCTGCGTGAAGAATTCA", start=2266659, alternate_bases=[
            "CACGCGAGTTGTAGATGATCGTTGAGTGGTCTTGCTTGGACTTCCATTTCATCTTTTCGACGCGCCAGGTCTCGCGGTCCTCCGGATCTGCGCCCGGTTTGAGTTGCACATCAAGGGGATACGGCTTGACCGACTCGTAGCCGACATGTAAGTCGGCTAGTTTCCGGCCGGCGCTGGCGAGCTGGTCGAAGCGTTCGCGGGTCTCCGGTGTTGGGATGTGCGGGAGCATCTTCTTGAGGTCAGCGGCGTATTTTGTGCGGTAGGCGGGGTCATGCAGCAGGCCGTAGACGTAGTAGAAGATGTCGTCTTTGGTGACTTGGTCGCCGATCGTGTCGCGGTAGAGCTTGAGGATGACGCCGGTGATGTTGTCGACGCGGCGGTAGCCGTGGTCGTCTACTTCGGCGTTGGTGGTGGACTCGAAATCGAGTTCGCCGTCACGTGGTTCGGTCTTCTCGTAGGTCCAGCGCGGGAAGAATTGACCGTTGCTTGAGCCCCAGAATGCGAGATCGGGGATAGCGTTTAGCATCAGACACGAGAAGGGCTTGTCTGAGCCCATGCCAACCACGTAGTAACCGACATTCCCGTGCTCCGGCGTCGGAAACATCGACGGAAGCTGGTAGGTACAGTTGTTGAGCTGCTGGTTGGGGTCGAGGTAGGCGTGCTCTTTCGTAAATGGTCGGTACGTGCCGAGCCGCATTCCCGCGGGAGCGAATTCGATGCGAATGCCTTGTGCCACTTGCCGCTTGTTGATGCGGTCCCAGCTGAACTTGGCCGAGTCCACGGTAATGAGGGCGTCAACCGGCGGGGTCTTGGCGTCCCTTCCGCGGATCTCGTTGATCCGGTCGACCTCCGAGTTGTAGAAGTCGATCGTGCGTCCGATGTTGGCCTCGAGCGCACCACGTGAAAAGTTGTAACACCACGCATCCCGGCTGGTCTTCAAGCCCGCGGAATAGTTCGCGAAGACACGTGTCACGTCAAGAGCAGCCTTCTTGTCGCCGATAACCGGCCACGCGCTGAACGCGTCGTCGCGTTGGTTGACCCAGTCACCGTGCAAGTTGGGTGTGACTGTCTGCCATTCCACCGTGTCGAGGTAGCCGTCGCCGACGATCCGCAACTTCTCCTCGCGACTCAGGTAATCGCCGATGTCGCGGTAAAGGACATCGCATGGCCCGCTGTGCTTCGGATCCTTGATGCCAAGGAAGATCGCCACCGTGTTGCGACTCCCCCCGCCAAAGACCTTGCCGCCTTCCTGGCGTGAGAGTTCCCCAGCTGTGCGCTGGTTCCCCCGCAGGTTGTACACATATACCGCCGCGTAGTCGTCGGCGAGCGACAACCGCATGCCGTCTGCCGTGTTGCCGTCTATGTACCCACCATTGGAGACGAATCCGACAACACCGTTGTCACCAATGCGGTCGGTCGCCCACCGGAACGCGCGAATATACGAGTCGTACAGGCTGTTCTTCAGCTGCGCCGTCGACCGCTTCGCGTACGTCTGCTCAATCCGCCCGTCCAACGTCGGATACTTCACGTTGGCGTTCAGGTCGTTCGCGCTGCTCTGCCCCACCGAGTACGGCGGATTCCCGATGATCACGCTGATCGGCGTCGCCAGCTGTCGCAAGATCCGAGCGTTGTTGTACGGGAACATGATCGCGTCCATCGAGTCCCCGGCTTCGGAAATCTGGAACGTGTCGGCCAGCGCCATCCCGGGGAACGGCTCATAGGCGTCGGCGTCGGCGGTCTTGCCCGCCAAAGCATGGTAGGTCGACTCGATGTTCACCGCGGCGATGTAGTACGCCAGCAGCATGATCTCGTTGGCGTGCAGCTCTTGCGAGTACTTTCGGGTGAGGTCGGCGGCCGTGATCAGGTCGGACTGCAGCAGCCGGGTAATGAATGTGCCCGTCCCGGCGAAGCCGTCCAGAATATGCACGCCCTCGTCGGTCAGCCCGCGCCCGAAATGCTTGCGCGACACGAAATCAGCCGCCCGCACAATGAAGTCCACGACCTCGACCGGCGTGTACACGATCCCCAGCGCCTCGGCCTGCTTCTTGAAGCCGATGCGGAAGAACTTCTCGTACAGCTCGGCGATCACCTGCTGCTTGCCCTCGGCGCTGGTGACCTCGCCGGCGCGCCGTCGCACCGATTCGTAAAAGCCTTCCAACCGAGCGGTTTCGGCCTCCAGGCCGGCACCCCCGACGGTGTCGACCATCTTCTGCATGGCCCGCGACACCGGGTTGTGCGACGCGAAGTCATGCCCGGCGAACAGCGCGTCGAACACCGGCTTGGTGATCAGGTGCTGCGAGAGCATGCTGATCGCGTCATCGGGGGTGATCGAGTCATTGAGGTTATCGCGCAGCCCGGCCAGGAACTGCTCGAACGCCGCCGCCGCCGTAGCGTCGGCGCCGCCGAGCAGGGCGTGGATACGGGTGGTCAGCGTCGCGGCGATGTCGGCGACATCGGCGGCCCACTGCTCCCAATAGGTCCGGGTGCCAACCTTGTCGACGATGCGCGCGTAGATCGCTTCCTGCCACTGCGACAACGAGAACATCGCCAACTGCTCCGCGACGGCGGGTCCCGCCTCGTCGGAGGTCGGCCCGATGTGACCGCCCAACAGCTTGTCGCTGCCTTCACCGGTCTTCGTCGGCTTCACGTTCAGCGCAATGCTGTTCACCATCGCGTCGAAGCGCTCGTCGTGCGACCGCAACGCGTTGAGGACCTGCCACACCACCTTGAACCGTTTGTTGTCGGCCAACGCGGCAGACGGCTCGACACCCTCGGGCACCGCCACCGGCAAGATGACGTACCCGTAGTCCTTGCCGGGCGACTTGCGCATCACCCGACCGACCGACTGCACCACGTCGACGATGGAATTGCGCGGATTCAGGAACAGCACCGCGTCCAGCGCGGGCACGTCGACCCCTTCGGAGAGGCAGCGGGCGTTGGACAGGATGCGGCATTCATCCTCGGCGACCACGCCTTTGAGCCAGGCCAGCTGTTCGTTGCGGACCAGCGCGTTGAACGTCCCGTCCACGTGGCGCACCG"])
        panel = self.pg.create(v)
        assert panel.ref == "TGGTGACGCGGGAGTAGAACGATCGCCAAGTGGTCGGTCTTGGCTGCCCACTTCATCCCCGGCGCCACCGGCAGGTCTCGCGGTCATCTCGACCAACGGAGGGCCGTCGGTGGTTCGTATCCGGCCAAGAACGGCGAGAACGGTTTGTGCCTCTATGCCAGGGTGAATGTCTCATCTCCCAGGCGGACGGTGATATCCAGTTCTCCGCCAAGAGCGGACACGTATTTGCGCAGTGTGTTGACCTGTGCGGAGCCGATGTCGCCGTTCTCGATGCTGGATACCCGGCTCTGCCGGATGTGCGCCAGCGCAGCCACCTGGACCTGGGTGAGTGACTGAGCCGCGCGCAGCTCCCGGAGCCGGAATGCCCGCACTTCATCGCGCATTCGTGCCTTGTGCCGGTCCACCGCCTCCCGGTTAACGGGACGTACGGCGTCCATGTCCCGTAGTGTCATCGCCATCGTGCCACTTACCCTTTCTTGCGCTTGCGCCTCTTTGGCTTCGTGTCCTCGAACTGTGCGAGATGTTCGGCAAACATCTCATCGGCCGCTTTGATCTTCTCGTCGTACCACTGGGTCCACCGCCCGGCCTTGTTACCGGCGGCCAGCATGATCGCCTGCCGCGCCGGGTCGAAGGCGAACAGAATGCGGACCTCGGACCGCCCTTGTGATCCTGGACGCAGCTCCTTCATGTTCTTGTGGCGCGACCCACGCACCGTGTCCACCAGAGGACAGCCAAGTGCGGGGCCCTCTTCCTCGAGAACCTCGATAGCTGCGAACACCAATTCGTAGGTCTCTCGGTCCAAGCCGTTGAGCCAGGCGGAGATGCGCTCCACATCCGCCGTCCACCCCACAGAGTCGCAGAGTAGCGCGATACGCGATATCACACAAGGGTGATATTCCTCCGGGTAAGAGCAGCGGGCGACGGGGCTACCGTCGAGGAAATGCCGGCAGGCGAGGACGGACTCTGCGCACCCGGGCCGTTGAAACAGTAGCCTGTGCCAGGCCGAGAATTCATCCCCACGTATGAGGCAGTACAGTGCGCCGCCGTGCGCGTTCTCCCATGGAACGTTCACGGGCTCCCGTGGATGACAGGCGTTTCATGAACGCCAGCGCCGCCGCAACCCGACCGAAAGCGGTTGACCCCAAGGAGAGCTGGAAGTCGAGGCCACCACCTTCGCCGCGGAGTTGCTCATGCCCGAGAGCGAGACTCGTCCCGAAATACGCCGGCTCGATTTCGGCAAGTTGCTCGAACTGAAGCGGGAATGGGCGTCGACCCGCTCGACCAGCCCCAGCCGGGTGACCAGCCCCAGCCGGGTGACCAGCCGATGCACCGCGGCGATCCCACCGAAGCCGGTGGCATCGATGTTGGCGCCGACCTCGTAGCGCACCGCGCCCGAACCCAGCATCGGCCTGGGCTGCGCCGCCCAGCGTCCAGCCCGCGCGTGCCGCGCCGCCACCCTGCGCCCTCGGCGTGTGATGTTTCGCCGACTCTGTTCATGGGTTATCTTCTTCACCACAAAGGCCTTTCCTGCTGGGCTGTGTTGAGGTCGCAAACCCAGCCAGGGTAAGGCCTTTGGCCTCTCCTACCCGGCCGACACGCTTACTGAAGGCCTAGTCTAGGCAGGCCATTCAATCTGCGGAATCGAAAAATTCGGTTCCAGCCTGCTCGTTTCCTTTCCGACAGCGATCTGACGTTGCGTAACGTCATTTGTACGGACTCTTTTAGCGGCATTGATTTCAGATGCCAACGCCGTCTGTGCTGTAGCGCCGATTGGCCGAAACTGTAAATTTGTATGATTATTTAAATCTTTGACGAACACGCGCCACAAACGTACTATCTCTTTGGCAAAGTCCACCGGCATCTCATTCAACGGTTTTGTTTGCGCGTGGTCGTCATATGTTGGTAACTGTGTAACCGGCCGCCTATCTTGCGCGTGCATCATATGACTATGAATCGGCCTTCTCCAGTGAAATTGATACAAGATCGATCCGATAAGCGGTACCTTGTACACAGTGCAATTGTAGTAATTCGCGTTTTGTCCTACGCTTGTATTCTGCGTGAAGAATTCAAACACG"
        assert panel.alts == [
            "GACCGCCGAGTGCGGCTGGATTGGATTTCACAAGGATGCCAATATCCGGCGCAACGCCGTCGAGCGACGGACGGTGCTCGACACGGGAGCCCGGCTATTCTGTGTGCCGCGGGCCGACATCCTGGCAGAGCAAGTCGCGGCACGGTATATTGCGTCCCTTGCGGCGATTGCCCGTGCCGCACGATTTCCGGGACCATTCATCTACACGGTTCACCCGAGCAAGATCGTTCGCGTGCTCTAGTCGTTCATCGCTCCGTTAACCGCCGGCGAGGCCGTCGACGATCTTCATGGTCTCGACGCTGACGGTGGTCACCTTCTTGATGAGGTCGACGATGTAGGTGGGATCGTCGTGTTCGTCGCACCAGTCGTTGGGGTCGTTGACGATGCCCGACGCTTTGTCGGTGGTGACGCGGTAGCGCTCGATGATCCAGCCGAGCGCCGAGCGGGAGCGAGCAGGTAGCGCTCGGCCTCGTCGGGAATGCCGGCGATGGTGACACGCGAGTTGTAGATGATCGTTGAGTGGTCTTGCTTGGACTTCCATTTCATCTTTTCGACGCGCCAGGTCTCGCGGTCCTCCGGATCTGCGCCCGGTTTGAGTTGCACATCAAGGGGATACGGCTTGACCGACTCGTAGCCGACATGTAAGTCGGCTAGTTTCCGGCCGGCGCTGGCGAGCTGGTCGAAGCGTTCGCGGGTCTCCGGTGTTGGGATGTGCGGGAGCATCTTCTTGAGGTCAGCGGCGTATTTTGTGCGGTAGGCGGGGTCATGCAGCAGGCCGTAGACGTAGTAGAAGATGTCGTCTTTGGTGACTTGGTCGCCGATCGTGTCGCGGTAGAGCTTGAGGATGACGCCGGTGATGTTGTCGACGCGGCGGTAGCCGTGGTCGTCTACTTCGGCGTTGGTGGTGGACTCGAAATCGAGTTCGCCGTCACGTGGTTCGGTCTTCTCGTAGGTCCAGCGCGGGAAGAATTGACCGTTGCTTGAGCCCCAGAATGCGAGATCGGGGATAGCGTTTAGCATCAGACACGAGAAGGGCTTGTCTGAGCCCATGCCAACCACGTAGTAACCGACATTCCCGTGCTCCGGCGTCGGAAACATCGACGGAAGCTGGTAGGTACAGTTGTTGAGCTGCTGGTTGGGGTCGAGGTAGGCGTGCTCTTTCGTAAATGGTCGGTACGTGCCGAGCCGCATTCCCGCGGGAGCGAATTCGATGCGAATGCCTTGTGCCACTTGCCGCTTGTTGATGCGGTCCCAGCTGAACTTGGCCGAGTCCACGGTAATGAGGGCGTCAACCGGCGGGGTCTTGGCGTCCCTTCCGCGGATCTCGTTGATCCGGTCGACCTCCGAGTTGTAGAAGTCGATCGTGCGTCCGATGTTGGCCTCGAGCGCACCACGTGAAAAGTTGTAACACCACGCATCCCGGCTGGTCTTCAAGCCCGCGGAATAGTTCGCGAAGACACGTGTCACGTCAAGAGCAGCCTTCTTGTCGCCGATAACCGGCCACGCGCTGAACGCGTCGTCGCGTTGGTTGACCCAGTCACCGTGCAAGTTGGGTGTGACTGTCTGCCATTCCACCGTGTCGAGGTAGCCGTCGCCGACGATCCGCAACTTCTCCTCGCGACTCAGGTAATCGCCGATGTCGCGGTAAAGGACATCGCATGGCCCGCTGTGCTTCGGATCCTTGATGCCAAGGAAGATCGCCACCGTGTTGCGACTCCCCCCGCCAAAGACCTTGCCGCCTTCCTGGCGTGAGAGTTCCCCAGCTGTGCGCTGGTTCCCCCGCAGGTTGTACACATATACCGCCGCGTAGTCGTCGGCGAGCGACAACCGCATGCCGTCTGCCGTGTTGCCGTCTATGTACCCACCATTGGAGACGAATCCGACAACACCGTTGTCACCAATGCGGTCGGTCGCCCACCGGAACGCGCGAATATACGAGTCGTACAGGCTGTTCTTCAGCTGCGCCGTCGACCGCTTCGCGTACGTCTGCTCAATCCGCCCGTCCAACGTCGGATACTTCACGTTGGCGTTCAGGTCGTTCGCGCTGCTCTGCCCCACCGAGTACGGCGGATTCCCGATGATCACGCTGATCGGCGTCGCCAGCTGTCGCAAGATCCGAGCGTTGTTGTACGGGAACATGATCGCGTCCATCGAGTCCCCGGCTTCGGAAATCTGGAACGTGTCGGCCAGCGCCATCCCGGGGAACGGCTCATAGGCGTCGGCGTCGGCGGTCTTGCCCGCCAAAGCATGGTAGGTCGACTCGATGTTCACCGCGGCGATGTAGTACGCCAGCAGCATGATCTCGTTGGCGTGCAGCTCTTGCGAGTACTTTCGGGTGAGGTCGGCGGCCGTGATCAGGTCGGACTGCAGCAGCCGGGTAATGAATGTGCCCGTCCCGGCGAAGCCGTCCAGAATATGCACGCCCTCGTCGGTCAGCCCGCGCCCGAAATGCTTGCGCGACACGAAATCAGCCGCCCGCACAATGAAGTCCACGACCTCGACCGGCGTGTACACGATCCCCAGCGCCTCGGCCTGCTTCTTGAAGCCGATGCGGAAGAACTTCTCGTACAGCTCGGCGATCACCTGCTGCTTGCCCTCGGCGCTGGTGACCTCGCCGGCGCGCCGTCGCACCGATTCGTAAAAGCCTTCCAACCGAGCGGTTTCGGCCTCCAGGCCGGCACCCCCGACGGTGTCGACCATCTTCTGCATGGCCCGCGACACCGGGTTGTGCGACGCGAAGTCATGCCCGGCGAACAGCGCGTCGAACACCGGCTTGGTGATCAGGTGCTGCGAGAGCATGCTGATCGCGTCATCGGGGGTGATCGAGTCATTGAGGTTATCGCGCAGCCCGGCCAGGAACTGCTCGAACGCCGCCGCCGCCGTAGCGTCGGCGCCGCCGAGCAGGGCGTGGATACGGGTGGTCAGCGTCGCGGCGATGTCGGCGACATCGGCGGCCCACTGCTCCCAATAGGTCCGGGTGCCAACCTTGTCGACGATGCGCGCGTAGATCGCTTCCTGCCACTGCGACAACGAGAACATCGCCAACTGCTCCGCGACGGCGGGTCCCGCCTCGTCGGAGGTCGGCCCGATGTGACCGCCCAACAGCTTGTCGCTGCCTTCACCGGTCTTCGTCGGCTTCACGTTCAGCGCAATGCTGTTCACCATCGCGTCGAAGCGCTCGTCGTGCGACCGCAACGCGTTGAGGACCTGCCACACCACCTTGAACCGTTTGTTGTCGGCCAACGCGGCAGACGGCTCGACACCCTCGGGCACCGCCACCGGCAAGATGACGTACCCGTAGTCCTTGCCGGGCGACTTGCGCATCACCCGACCGACCGACTGCACCACGTCGACGATGGAATTGCGCGGATTCAGGAACAGCACCGCGTCCAGCGCGGGCACGTCGACCCCTTCGGAGAGGCAGCGGGCGTTGGACAGGATGCGGCATTCATCCTCGGCGACCACGCCTTTGAGCCAGGCCAGCTGTTCGTTGCGGACCAGCGCGTTGAACGTCCCGTCCACGTGGCGCACCGAACACGCCAGGCCCGGGCCGTCGTCAACCAATTCGCGGTATGCCTCAACCACTTTCGGGAACAGCTCGGCAACCTGCTTGGACGTCTTGATGTCCTTGGCGAACGCCACCGCCCGACGCATCGGCGGCTCACCGGCGACAATGCCGGTACCGGACCGCTTGGCCAGGCCATTCCAGCAGCCGACGATCTTGGAGGCGTCGTCGAGCATCAGCTCGCCGGAAACCCCGGAGAGTTCCTGCTGCAACCGGGGCGCGATCACGCCCTGATCGACGGTGAGCACCATCACCTTGTAGTCGGTGAGCAGCCCGCGCTCCACCGCCTCGCCGAACGACAGCCGGTGAAACTCCGGCCCGAACGTCAGCTCGTCGTCCATCGACACCAACTCGGCGGAGTGCTGGTCGGCCCTGTCCTTGATGCTCTCGGTGAAAATCCTTGGCGTGGCGGTCATATACAGCCGCCGGGCCGCCTTCAGATACTGACCGTCGTGCACCCGC"]
