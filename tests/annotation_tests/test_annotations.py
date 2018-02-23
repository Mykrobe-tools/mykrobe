import copy

from mykrobe.annotation.genes import Region
from mykrobe.annotation.genes import Gene
from mykrobe.annotation.genes import GeneAminoAcidChangeToDNAVariants

from mykrobe.panelgeneration import AlleleGenerator
from ga4ghmongo.schema import Variant
from ga4ghmongo.schema import VariantSet
from ga4ghmongo.schema import Reference
from ga4ghmongo.schema import ReferenceSet
from mykrobe.utils import split_var_name

from Bio import SeqIO
from Bio.Seq import Seq
from mongoengine import connect
DB = connect('atlas-test')
import pytest


class TestRegions():

    def setup(self):
        DB.drop_database('atlas-test')
        with open("src/mykrobe/data/NC_000962.3.fasta", 'r') as infile:
            self.reference_seq = list(SeqIO.parse(infile, "fasta"))[0].seq

        self.gm = GeneAminoAcidChangeToDNAVariants(
            reference="src/mykrobe/data/NC_000962.3.fasta",
            genbank="src/mykrobe/data/NC_000962.3.gb")
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.variant_set = VariantSet.create_and_save(
            name="this_vcf_file",
            reference_set=self.reference_set)
        self.variant_sets = [self.variant_set]
        self.reference_id = Reference().create_and_save(
            name="ref",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])

    def test_simple_gene(self):
        g = Gene(
            name="rpoB",
            reference=self.reference_seq,
            start=759807,
            end=763325)
        assert g.name == "rpoB"
        assert g.forward
        assert g.strand == "forward"
        assert g.seq == "TTGGCAGATTCCCGCCAGAGCAAAACAGCCGCTAGTCCTAGTCCGAGTCGCCCGCAAAGTTCCTCGAATAACTCCGTACCCGGAGCGCCAAACCGGGTCTCCTTCGCTAAGCTGCGCGAACCACTTGAGGTTCCGGGACTCCTTGACGTCCAGACCGATTCGTTCGAGTGGCTGATCGGTTCGCCGCGCTGGCGCGAATCCGCCGCCGAGCGGGGTGATGTCAACCCAGTGGGTGGCCTGGAAGAGGTGCTCTACGAGCTGTCTCCGATCGAGGACTTCTCCGGGTCGATGTCGTTGTCGTTCTCTGACCCTCGTTTCGACGATGTCAAGGCACCCGTCGACGAGTGCAAAGACAAGGACATGACGTACGCGGCTCCACTGTTCGTCACCGCCGAGTTCATCAACAACAACACCGGTGAGATCAAGAGTCAGACGGTGTTCATGGGTGACTTCCCGATGATGACCGAGAAGGGCACGTTCATCATCAACGGGACCGAGCGTGTGGTGGTCAGCCAGCTGGTGCGGTCGCCCGGGGTGTACTTCGACGAGACCATTGACAAGTCCACCGACAAGACGCTGCACAGCGTCAAGGTGATCCCGAGCCGCGGCGCGTGGCTCGAGTTTGACGTCGACAAGCGCGACACCGTCGGCGTGCGCATCGACCGCAAACGCCGGCAACCGGTCACCGTGCTGCTCAAGGCGCTGGGCTGGACCAGCGAGCAGATTGTCGAGCGGTTCGGGTTCTCCGAGATCATGCGATCGACGCTGGAGAAGGACAACACCGTCGGCACCGACGAGGCGCTGTTGGACATCTACCGCAAGCTGCGTCCGGGCGAGCCCCCGACCAAAGAGTCAGCGCAGACGCTGTTGGAAAACTTGTTCTTCAAGGAGAAGCGCTACGACCTGGCCCGCGTCGGTCGCTATAAGGTCAACAAGAAGCTCGGGCTGCATGTCGGCGAGCCCATCACGTCGTCGACGCTGACCGAAGAAGACGTCGTGGCCACCATCGAATATCTGGTCCGCTTGCACGAGGGTCAGACCACGATGACCGTTCCGGGCGGCGTCGAGGTGCCGGTGGAAACCGACGACATCGACCACTTCGGCAACCGCCGCCTGCGTACGGTCGGCGAGCTGATCCAAAACCAGATCCGGGTCGGCATGTCGCGGATGGAGCGGGTGGTCCGGGAGCGGATGACCACCCAGGACGTGGAGGCGATCACACCGCAGACGTTGATCAACATCCGGCCGGTGGTCGCCGCGATCAAGGAGTTCTTCGGCACCAGCCAGCTGAGCCAATTCATGGACCAGAACAACCCGCTGTCGGGGTTGACCCACAAGCGCCGACTGTCGGCGCTGGGGCCCGGCGGTCTGTCACGTGAGCGTGCCGGGCTGGAGGTCCGCGACGTGCACCCGTCGCACTACGGCCGGATGTGCCCGATCGAAACCCCTGAGGGGCCCAACATCGGTCTGATCGGCTCGCTGTCGGTGTACGCGCGGGTCAACCCGTTCGGGTTCATCGAAACGCCGTACCGCAAGGTGGTCGACGGCGTGGTTAGCGACGAGATCGTGTACCTGACCGCCGACGAGGAGGACCGCCACGTGGTGGCACAGGCCAATTCGCCGATCGATGCGGACGGTCGCTTCGTCGAGCCGCGCGTGCTGGTCCGCCGCAAGGCGGGCGAGGTGGAGTACGTGCCCTCGTCTGAGGTGGACTACATGGACGTCTCGCCCCGCCAGATGGTGTCGGTGGCCACCGCGATGATTCCCTTCCTGGAGCACGACGACGCCAACCGTGCCCTCATGGGGGCAAACATGCAGCGCCAGGCGGTGCCGCTGGTCCGTAGCGAGGCCCCGCTGGTGGGCACCGGGATGGAGCTGCGCGCGGCGATCGACGCCGGCGACGTCGTCGTCGCCGAAGAAAGCGGCGTCATCGAGGAGGTGTCGGCCGACTACATCACTGTGATGCACGACAACGGCACCCGGCGTACCTACCGGATGCGCAAGTTTGCCCGGTCCAACCACGGCACTTGCGCCAACCAGTGCCCCATCGTGGACGCGGGCGACCGAGTCGAGGCCGGTCAGGTGATCGCCGACGGTCCCTGTACTGACGACGGCGAGATGGCGCTGGGCAAGAACCTGCTGGTGGCCATCATGCCGTGGGAGGGCCACAACTACGAGGACGCGATCATCCTGTCCAACCGCCTGGTCGAAGAGGACGTGCTCACCTCGATCCACATCGAGGAGCATGAGATCGATGCTCGCGACACCAAGCTGGGTGCGGAGGAGATCACCCGCGACATCCCGAACATCTCCGACGAGGTGCTCGCCGACCTGGATGAGCGGGGCATCGTGCGCATCGGTGCCGAGGTTCGCGACGGGGACATCCTGGTCGGCAAGGTCACCCCGAAGGGTGAGACCGAGCTGACGCCGGAGGAGCGGCTGCTGCGTGCCATCTTCGGTGAGAAGGCCCGCGAGGTGCGCGACACTTCGCTGAAGGTGCCGCACGGCGAATCCGGCAAGGTGATCGGCATTCGGGTGTTTTCCCGCGAGGACGAGGACGAGTTGCCGGCCGGTGTCAACGAGCTGGTGCGTGTGTATGTGGCTCAGAAACGCAAGATCTCCGACGGTGACAAGCTGGCCGGCCGGCACGGCAACAAGGGCGTGATCGGCAAGATCCTGCCGGTTGAGGACATGCCGTTCCTTGCCGACGGCACCCCGGTGGACATTATTTTGAACACCCACGGCGTGCCGCGACGGATGAACATCGGCCAGATTTTGGAGACCCACCTGGGTTGGTGTGCCCACAGCGGCTGGAAGGTCGACGCCGCCAAGGGGGTTCCGGACTGGGCCGCCAGGCTGCCCGACGAACTGCTCGAGGCGCAGCCGAACGCCATTGTGTCGACGCCGGTGTTCGACGGCGCCCAGGAGGCCGAGCTGCAGGGCCTGTTGTCGTGCACGCTGCCCAACCGCGACGGTGACGTGCTGGTCGACGCCGACGGCAAGGCCATGCTCTTCGACGGGCGCAGCGGCGAGCCGTTCCCGTACCCGGTCACGGTTGGCTACATGTACATCATGAAGCTGCACCACCTGGTGGACGACAAGATCCACGCCCGCTCCACCGGGCCGTACTCGATGATCACCCAGCAGCCGCTGGGCGGTAAGGCGCAGTTCGGTGGCCAGCGGTTCGGGGAGATGGAGTGCTGGGCCATGCAGGCCTACGGTGCTGCCTACACCCTGCAGGAGCTGTTGACCATCAAGTCCGATGACACCGTCGGCCGCGTCAAGGTGTACGAGGCGATCGTCAAGGGTGAGAACATCCCGGAGCCGGGCATCCCCGAGTCGTTCAAGGTGCTGCTCAAAGAACTGCAGTCGCTGTGCCTCAACGTCGAGGTGCTATCGAGTGACGGTGCGGCGATCGAACTGCGCGAAGGTGAGGACGAGGACCTGGAGCGGGCCGCGGCCAACCTGGGAATCAATCTGTCCCGCAACGAATCCGCAAGTGTCGAGGATCTTGCGTAA"
        assert g.prot == "LADSRQSKTAASPSPSRPQSSSNNSVPGAPNRVSFAKLREPLEVPGLLDVQTDSFEWLIGSPRWRESAAERGDVNPVGGLEEVLYELSPIEDFSGSMSLSFSDPRFDDVKAPVDECKDKDMTYAAPLFVTAEFINNNTGEIKSQTVFMGDFPMMTEKGTFIINGTERVVVSQLVRSPGVYFDETIDKSTDKTLHSVKVIPSRGAWLEFDVDKRDTVGVRIDRKRRQPVTVLLKALGWTSEQIVERFGFSEIMRSTLEKDNTVGTDEALLDIYRKLRPGEPPTKESAQTLLENLFFKEKRYDLARVGRYKVNKKLGLHVGEPITSSTLTEEDVVATIEYLVRLHEGQTTMTVPGGVEVPVETDDIDHFGNRRLRTVGELIQNQIRVGMSRMERVVRERMTTQDVEAITPQTLINIRPVVAAIKEFFGTSQLSQFMDQNNPLSGLTHKRRLSALGPGGLSRERAGLEVRDVHPSHYGRMCPIETPEGPNIGLIGSLSVYARVNPFGFIETPYRKVVDGVVSDEIVYLTADEEDRHVVAQANSPIDADGRFVEPRVLVRRKAGEVEYVPSSEVDYMDVSPRQMVSVATAMIPFLEHDDANRALMGANMQRQAVPLVRSEAPLVGTGMELRAAIDAGDVVVAEESGVIEEVSADYITVMHDNGTRRTYRMRKFARSNHGTCANQCPIVDAGDRVEAGQVIADGPCTDDGEMALGKNLLVAIMPWEGHNYEDAIILSNRLVEEDVLTSIHIEEHEIDARDTKLGAEEITRDIPNISDEVLADLDERGIVRIGAEVRDGDILVGKVTPKGETELTPEERLLRAIFGEKAREVRDTSLKVPHGESGKVIGIRVFSREDEDELPAGVNELVRVYVAQKRKISDGDKLAGRHGNKGVIGKILPVEDMPFLADGTPVDIILNTHGVPRRMNIGQILETHLGWCAHSGWKVDAAKGVPDWAARLPDELLEAQPNAIVSTPVFDGAQEAELQGLLSCTLPNRDGDVLVDADGKAMLFDGRSGEPFPYPVTVGYMYIMKLHHLVDDKIHARSTGPYSMITQQPLGGKAQFGGQRFGEMECWAMQAYGAAYTLQELLTIKSDDTVGRVKVYEAIVKGENIPEPGIPESFKVLLKELQSLCLNVEVLSSDGAAIELREGEDEDLERAAANLGINLSRNESASVEDLA"

    def test_reverse_gene(self):
        g = Gene(
            name="gidB",
            reference=self.reference_seq,
            start=4407528,
            end=4408202,
            forward=False)
        assert g.name == "gidB"
        assert g.forward is False
        assert g.strand == "reverse"
        assert g.seq == "ATGTCTCCGATCGAGCCCGCGGCGTCTGCGATCTTCGGACCGCGGCTTGGCCTTGCTCGGCGGTACGCCGAAGCGTTGGCGGGACCCGGTGTGGAGCGGGGGCTGGTGGGACCCCGCGAAGTCGGTAGGCTATGGGACCGGCATCTACTGAACTGCGCCGTGATCGGTGAGCTCCTCGAACGCGGTGACCGGGTCGTGGATATCGGTAGCGGAGCCGGGTTGCCGGGCGTGCCATTGGCGATAGCGCGGCCGGACCTCCAGGTAGTTCTCCTAGAACCGCTACTGCGCCGCACCGAGTTTCTTCGAGAGATGGTGACAGATCTGGGCGTGGCCGTTGAGATCGTGCGGGGGCGCGCCGAGGAGTCCTGGGTGCAGGACCAATTGGGCGGCAGCGACGCTGCGGTGTCACGGGCGGTGGCCGCGTTGGACAAGTTGACGAAATGGAGCATGCCGTTGATACGGCCGAACGGGCGAATGCTCGCCATCAAAGGCGAGCGGGCTCACGACGAAGTACGGGAGCACCGGCGTGTGATGATCGCATCGGGCGCGGTTGATGTCAGGGTGGTGACATGTGGCGCGAACTATTTGCGTCCGCCCGCGACCGTGGTGTTCGCACGACGTGGAAAGCAGATCGCCCGAGGGTCGGCACGGATGGCGAGTGGAGGGACGGCGTGA"
        assert g.prot == "MSPIEPAASAIFGPRLGLARRYAEALAGPGVERGLVGPREVGRLWDRHLLNCAVIGELLERGDRVVDIGSGAGLPGVPLAIARPDLQVVLLEPLLRRTEFLREMVTDLGVAVEIVRGRAEESWVQDQLGGSDAAVSRAVAALDKLTKWSMPLIRPNGRMLAIKGERAHDEVREHRRVMIASGAVDVRVVTCGANYLRPPATVVFARRGKQIARGSARMASGGTA"

    def test_reverse_gene2(self):
        g = Gene(
            name="katG",
            reference=self.reference_seq,
            start=2153889,
            end=2156111,
            forward=False)
        assert g.name == "katG"
        assert g.forward is False
        assert g.strand == "reverse"
        assert g.seq == "GTGCCCGAGCAACACCCACCCATTACAGAAACCACCACCGGAGCCGCTAGCAACGGCTGTCCCGTCGTGGGTCATATGAAATACCCCGTCGAGGGCGGCGGAAACCAGGACTGGTGGCCCAACCGGCTCAATCTGAAGGTACTGCACCAAAACCCGGCCGTCGCTGACCCGATGGGTGCGGCGTTCGACTATGCCGCGGAGGTCGCGACCATCGACGTTGACGCCCTGACGCGGGACATCGAGGAAGTGATGACCACCTCGCAGCCGTGGTGGCCCGCCGACTACGGCCACTACGGGCCGCTGTTTATCCGGATGGCGTGGCACGCTGCCGGCACCTACCGCATCCACGACGGCCGCGGCGGCGCCGGGGGCGGCATGCAGCGGTTCGCGCCGCTTAACAGCTGGCCCGACAACGCCAGCTTGGACAAGGCGCGCCGGCTGCTGTGGCCGGTCAAGAAGAAGTACGGCAAGAAGCTCTCATGGGCGGACCTGATTGTTTTCGCCGGCAACTGCGCGCTGGAATCGATGGGCTTCAAGACGTTCGGGTTCGGCTTCGGCCGGGTCGACCAGTGGGAGCCCGATGAGGTCTATTGGGGCAAGGAAGCCACCTGGCTCGGCGATGAGCGTTACAGCGGTAAGCGGGATCTGGAGAACCCGCTGGCCGCGGTGCAGATGGGGCTGATCTACGTGAACCCGGAGGGGCCGAACGGCAACCCGGACCCCATGGCCGCGGCGGTCGACATTCGCGAGACGTTTCGGCGCATGGCCATGAACGACGTCGAAACAGCGGCGCTGATCGTCGGCGGTCACACTTTCGGTAAGACCCATGGCGCCGGCCCGGCCGATCTGGTCGGCCCCGAACCCGAGGCTGCTCCGCTGGAGCAGATGGGCTTGGGCTGGAAGAGCTCGTATGGCACCGGAACCGGTAAGGACGCGATCACCAGCGGCATCGAGGTCGTATGGACGAACACCCCGACGAAATGGGACAACAGTTTCCTCGAGATCCTGTACGGCTACGAGTGGGAGCTGACGAAGAGCCCTGCTGGCGCTTGGCAATACACCGCCAAGGACGGCGCCGGTGCCGGCACCATCCCGGACCCGTTCGGCGGGCCAGGGCGCTCCCCGACGATGCTGGCCACTGACCTCTCGCTGCGGGTGGATCCGATCTATGAGCGGATCACGCGTCGCTGGCTGGAACACCCCGAGGAATTGGCCGACGAGTTCGCCAAGGCCTGGTACAAGCTGATCCACCGAGACATGGGTCCCGTTGCGAGATACCTTGGGCCGCTGGTCCCCAAGCAGACCCTGCTGTGGCAGGATCCGGTCCCTGCGGTCAGCCACGACCTCGTCGGCGAAGCCGAGATTGCCAGCCTTAAGAGCCAGATCCGGGCATCGGGATTGACTGTCTCACAGCTAGTTTCGACCGCATGGGCGGCGGCGTCGTCGTTCCGTGGTAGCGACAAGCGCGGCGGCGCCAACGGTGGTCGCATCCGCCTGCAGCCACAAGTCGGGTGGGAGGTCAACGACCCCGACGGGGATCTGCGCAAGGTCATTCGCACCCTGGAAGAGATCCAGGAGTCATTCAACTCCGCGGCGCCGGGGAACATCAAAGTGTCCTTCGCCGACCTCGTCGTGCTCGGTGGCTGTGCCGCCATAGAGAAAGCAGCAAAGGCGGCTGGCCACAACATCACGGTGCCCTTCACCCCGGGCCGCACGGATGCGTCGCAGGAACAAACCGACGTGGAATCCTTTGCCGTGCTGGAGCCCAAGGCAGATGGCTTCCGAAACTACCTCGGAAAGGGCAACCCGTTGCCGGCCGAGTACATGCTGCTCGACAAGGCGAACCTGCTTACGCTCAGTGCCCCTGAGATGACGGTGCTGGTAGGTGGCCTGCGCGTCCTCGGCGCAAACTACAAGCGCTTACCGCTGGGCGTGTTCACCGAGGCCTCCGAGTCACTGACCAACGACTTCTTCGTGAACCTGCTCGACATGGGTATCACCTGGGAGCCCTCGCCAGCAGATGACGGGACCTACCAGGGCAAGGATGGCAGTGGCAAGGTGAAGTGGACCGGCAGCCGCGTGGACCTGGTCTTCGGGTCCAACTCGGAGTTGCGGGCGCTTGTCGAGGTCTATGGCGCCGATGACGCGCAGCCGAAGTTCGTGCAGGACTTCGTCGCTGCCTGGGACAAGGTGATGAACCTCGACAGGTTCGACGTGCGCTGA"
        assert g.prot == "VPEQHPPITETTTGAASNGCPVVGHMKYPVEGGGNQDWWPNRLNLKVLHQNPAVADPMGAAFDYAAEVATIDVDALTRDIEEVMTTSQPWWPADYGHYGPLFIRMAWHAAGTYRIHDGRGGAGGGMQRFAPLNSWPDNASLDKARRLLWPVKKKYGKKLSWADLIVFAGNCALESMGFKTFGFGFGRVDQWEPDEVYWGKEATWLGDERYSGKRDLENPLAAVQMGLIYVNPEGPNGNPDPMAAAVDIRETFRRMAMNDVETAALIVGGHTFGKTHGAGPADLVGPEPEAAPLEQMGLGWKSSYGTGTGKDAITSGIEVVWTNTPTKWDNSFLEILYGYEWELTKSPAGAWQYTAKDGAGAGTIPDPFGGPGRSPTMLATDLSLRVDPIYERITRRWLEHPEELADEFAKAWYKLIHRDMGPVARYLGPLVPKQTLLWQDPVPAVSHDLVGEAEIASLKSQIRASGLTVSQLVSTAWAAASSFRGSDKRGGANGGRIRLQPQVGWEVNDPDGDLRKVIRTLEEIQESFNSAAPGNIKVSFADLVVLGGCAAIEKAAKAAGHNITVPFTPGRTDASQEQTDVESFAVLEPKADGFRNYLGKGNPLPAEYMLLDKANLLTLSAPEMTVLVGGLRVLGANYKRLPLGVFTEASESLTNDFFVNLLDMGITWEPSPADDGTYQGKDGSGKVKWTGSRVDLVFGSNSELRALVEVYGADDAQPKFVQDFVAAWDKVMNLDRFDVR"

    def test_get_codon(self):
        g = Gene(
            name="rpoB",
            reference=self.reference_seq,
            start=759807,
            end=763325)
        with pytest.raises(ValueError):
            g.get_codon(1173)
        assert g.get_codon(2) == "GCA"
        assert g.get_codon(3) == "GAT"
        assert g.get_reference_position(1) == 759807
        assert g.seq[0] == self.reference_seq[759806]
        assert g.get_reference_position(-1) == 759806

    def test_get_codon_reverse(self):
        g = Gene(
            name="gidB",
            reference=self.reference_seq,
            start=4407528,
            end=4408202,
            forward=False)
        with pytest.raises(ValueError):
            g.get_codon(225)
        assert g.get_codon(2) == "TCT"
        assert g.get_codon(3) == "CCG"
        assert g.get_reference_position(1) == 4408202
        assert g.get_reference_position(2) == 4408201
        assert g.get_reference_position(-1) == 4408203
        assert g.get_reference_position(-2) == 4408204

    def test_gene_muts(self):
        self.gm = GeneAminoAcidChangeToDNAVariants(
            reference="src/mykrobe/data/NC_000962.3.fasta",
            genbank="src/mykrobe/data/NC_000962.3.gb")
        assert self.gm.get_alts("K") == ['AAA', 'AAG']
        # GAT -> ['GCA', 'GCT', 'GCC', 'GCG'], positions 759813,14,15
        assert sorted(self.gm.get_variant_names("rpoB", "D3A")) == sorted(
            ['GAT759813GCA', 'GAT759813GCT', 'GAT759813GCC', 'GAT759813GCG'])
        # GAT -> ['GCA', 'GCT', 'GCC', 'GCG'], positions 759813,14,15
        assert sorted(self.gm.get_variant_names("rpoB",
                                                "D3X")) == sorted(['GAT759813GCA',
                                                                   'GAT759813GCT',
                                                                   'GAT759813GCC',
                                                                   'GAT759813GCG',
                                                                   'GAT759813TGT',
                                                                   'GAT759813TGC',
                                                                   'GAT759813GAA',
                                                                   'GAT759813GAG',
                                                                   'GAT759813GGA',
                                                                   'GAT759813GGT',
                                                                   'GAT759813GGC',
                                                                   'GAT759813GGG',
                                                                   'GAT759813TTT',
                                                                   'GAT759813TTC',
                                                                   'GAT759813ATA',
                                                                   'GAT759813ATT',
                                                                   'GAT759813ATC',
                                                                   'GAT759813CAT',
                                                                   'GAT759813CAC',
                                                                   'GAT759813AAA',
                                                                   'GAT759813AAG',
                                                                   'GAT759813ATG',
                                                                   'GAT759813TTA',
                                                                   'GAT759813TTG',
                                                                   'GAT759813CTA',
                                                                   'GAT759813CTT',
                                                                   'GAT759813CTC',
                                                                   'GAT759813CTG',
                                                                   'GAT759813AAT',
                                                                   'GAT759813AAC',
                                                                   'GAT759813CAA',
                                                                   'GAT759813CAG',
                                                                   'GAT759813CCA',
                                                                   'GAT759813CCT',
                                                                   'GAT759813CCC',
                                                                   'GAT759813CCG',
                                                                   'GAT759813AGT',
                                                                   'GAT759813AGC',
                                                                   'GAT759813TCA',
                                                                   'GAT759813TCT',
                                                                   'GAT759813TCC',
                                                                   'GAT759813TCG',
                                                                   'GAT759813AGA',
                                                                   'GAT759813AGG',
                                                                   'GAT759813CGA',
                                                                   'GAT759813CGT',
                                                                   'GAT759813CGC',
                                                                   'GAT759813CGG',
                                                                   'GAT759813ACA',
                                                                   'GAT759813ACT',
                                                                   'GAT759813ACC',
                                                                   'GAT759813ACG',
                                                                   'GAT759813TGG',
                                                                   'GAT759813GTA',
                                                                   'GAT759813GTT',
                                                                   'GAT759813GTC',
                                                                   'GAT759813GTG',
                                                                   'GAT759813TAT',
                                                                   'GAT759813TAC'])

    def test_gene_muts2(self):
        self.gm = GeneAminoAcidChangeToDNAVariants(
            reference="src/mykrobe/data/NC_000962.3.fasta",
            genbank="src/mykrobe/data/NC_000962.3.gb")
        assert self.gm.get_alts("K") == ['AAA', 'AAG']
        # AGC -> ['CTT', 'CTC', 'CTA', 'CTG']
    #   # GAG -> ['GCA', 'GCT', 'GCC', 'GCG']
        # RC : CTC -> ['TGC',...] position2156103
        assert sorted(self.gm.get_variant_names("katG", "E3A")) == sorted(
            ['CTC2156103TGC', 'CTC2156103AGC', 'CTC2156103GGC', 'CTC2156103CGC'])

    def test_make_variant_panel(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("rpoB")
        for var in self.gm.get_variant_names("rpoB", "D3A"):
            ref, start, alt = split_var_name(var)
            v = Variant.create(
                variant_sets=self.variant_sets,
                reference=self.reference_id,
                reference_bases=ref,
                start=start,
                alternate_bases=[alt])
            panel = ag.create(v)
            for alt in panel.alts:
                seq = copy.copy(str(gene.seq))
                seq = seq.replace(panel.ref[25:], alt[25:])
                assert seq != str(gene.seq)
                assert Seq(seq).translate()[2] == "A"

    def test_make_variant_panel2(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("katG")
        for var in self.gm.get_variant_names("katG", "E3A"):
            ref, start, alt = split_var_name(var)
            v = Variant.create(
                variant_sets=self.variant_sets,
                reference=self.reference_id,
                reference_bases=ref,
                start=start,
                alternate_bases=[alt])
            panel = ag.create(v)
            for alt in panel.alts:
                seq = copy.copy(str(gene.seq.reverse_complement()))
                seq = seq.replace(panel.ref[:39], alt[:39])
                assert seq != str(gene.seq)
                assert Seq(seq).reverse_complement().translate()[2] == "A"

    def test_make_variant_panel3(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("katG")
        for var in self.gm.get_variant_names("katG", "S315L"):
            ref, start, alt = split_var_name(var)
            v = Variant.create(
                variant_sets=self.variant_sets,
                reference=self.reference_id,
                reference_bases=ref,
                start=start,
                alternate_bases=[alt])
            panel = ag.create(v)
            for alt in panel.alts:
                seq = copy.copy(str(gene.seq.reverse_complement()))
                seq = seq.replace(panel.ref, alt)
                assert seq != str(gene.seq)
                assert Seq(seq).reverse_complement().translate()[314] == "L"

    def test_make_variant_panel4(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("katG")
        for var in self.gm.get_variant_names("katG", "W90R"):
            ref, start, alt = split_var_name(var)
            v = Variant.create(
                variant_sets=self.variant_sets,
                reference=self.reference_id,
                reference_bases=ref,
                start=start,
                alternate_bases=[alt])
            panel = ag.create(v)
            for alt in panel.alts:
                seq = copy.copy(str(gene.seq.reverse_complement()))
                seq = seq.replace(panel.ref, alt)
                assert seq != str(gene.seq)
                assert Seq(seq).reverse_complement().translate()[89] == "R"

    def test_make_variant_panel5(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("gyrA")
        for var in self.gm.get_variant_names("gyrA", "D94X"):
            ref, start, alt = split_var_name(var)
            v = Variant.create(
                variant_sets=self.variant_sets,
                reference=self.reference_id,
                reference_bases=ref,
                start=start,
                alternate_bases=[alt])
            panel = ag.create(v)
            for alt in panel.alts:
                seq = copy.copy(str(gene.seq))
                seq = seq.replace(panel.ref, alt)
                assert Seq(seq).translate()[93] != "D"

    def test_make_variant_panel6(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("pncA")
        variants = list(self.gm.get_variant_names(
            "pncA", "CAG28TAA", protein_coding_var=False))
        assert len(variants) == 1
        var = variants[0]
        ref, start, alt = split_var_name(var)
        assert ref == 'CTG'
        assert start == 2289212
        assert alt == 'TTA'
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt])
        panel = ag.create(v)
        assert len(panel.alts) == 1
        alt = panel.alts[0]
        # the panel ref/alt seqs go past the end of the gene,
        # so can't comparie against gene sequence. Need to get
        # subseq from the reference seq
        panel_ref_start = self.reference_seq.find(panel.ref)
        assert panel_ref_start < start < panel_ref_start + len(panel.ref)
        seq = str(
            self.reference_seq[panel_ref_start:panel_ref_start + len(panel.ref)])
        assert seq == panel.ref
        assert alt == seq[:30] + 'TTA' + seq[33:]

    def test_make_variant_panel7(self):
        # Test DNA change upstream of a gene on the reverse
        # strand. The variant G-10A is in "gene space", ie
        # 10 bases upstream of eis is the nucleotide G on the
        # reverse strand. That position is 2715342 in the genome,
        # and is C on the forwards strand.
        # Here's a diagram:
        #             | <- This C is at -10 in "gene space", so variant G-10A has ref=G
        #             |    ref coord is 2715342, and variant in "ref space" is C2715342T
        # CACAGAATCCGACTGTGGCATATGCCGC
        #   |
        #   | <- C = last nucleotide of gene, at 2715332
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("eis")
        variants = list(self.gm.get_variant_names(
            "eis", "G-10A", protein_coding_var=False))
        assert len(variants) == 1
        var = variants[0]
        ref, start, alt = split_var_name(var)
        assert ref == 'C'
        assert start == 2715342
        assert alt == 'T'
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt])
        panel = ag.create(v)
        assert len(panel.alts) == 1
        alt = panel.alts[0]
        # the panel ref/alt seqs go past the end of the gene,
        # so can't comparie against gene sequence. Need to get
        # subseq from the reference seq
        panel_ref_start = self.reference_seq.find(panel.ref)
        assert panel_ref_start < start < panel_ref_start + len(panel.ref)
        seq = str(
            self.reference_seq[panel_ref_start:panel_ref_start + len(panel.ref)])
        assert seq == panel.ref
        assert alt == seq[:30] + 'T' + seq[31:]

    def test_make_variant_panel8(self):
        ag = AlleleGenerator("src/mykrobe/data/NC_000962.3.fasta")
        gene = self.gm.get_gene("eis")
        variants = list(self.gm.get_variant_names(
            "eis", "TG-1T", protein_coding_var=False))
        assert len(variants) == 1
        var = variants[0]
        ref, start, alt = split_var_name(var)
        assert ref == 'CA'
        assert start == 2715332
        assert alt == 'A'
        v = Variant.create(
            variant_sets=self.variant_sets,
            reference=self.reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt])
        panel = ag.create(v)
        assert len(panel.alts) == 1
        alt = panel.alts[0]
        # the panel ref/alt seqs go past the end of the gene,
        # so can't comparie against gene sequence. Need to get
        # subseq from the reference seq
        panel_ref_start = self.reference_seq.find(panel.ref)
        assert panel_ref_start < start < panel_ref_start + len(panel.ref)
        seq = str(
            self.reference_seq[panel_ref_start:panel_ref_start + len(panel.ref)])
        assert seq == panel.ref
        assert alt[:-1] == seq[:30] + seq[31:]
