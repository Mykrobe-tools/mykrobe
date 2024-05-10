import copy
import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from mongoengine import connect

from mykrobe.annotation.genes import Gene
from mykrobe.annotation.genes import GeneAminoAcidChangeToDNAVariants
from mykrobe.probes import AlleleGenerator
from mykrobe.utils import split_var_name
from mykrobe.variants.schema.models import Reference
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantSet

DATA_DIR = os.path.join("tests", "ref_data")


@pytest.fixture(autouse=True)
def db_setup_teardown():
    DB = connect("mykrobe-test", uuidRepresentation="pythonLegacy")
    DB.drop_database("mykrobe-test")
    yield
    DB.drop_database("mykrobe-test")


@pytest.fixture()
def reference_seq():
    with open(f"{DATA_DIR}/NC_000962.3.fasta", "r") as infile:
        reference_seq = list(SeqIO.parse(infile, "fasta"))[0].seq
    return reference_seq


@pytest.fixture()
def reference_id_and_variant_sets(reference_seq):
    reference_set = ReferenceSet().create_and_save(name="ref_set")
    variant_set = VariantSet.create_and_save(
        name="this_vcf_file", reference_set=reference_set
    )
    variant_sets = [variant_set]
    reference_id = Reference().create_and_save(
        name="ref", md5checksum="sre", reference_sets=[reference_set]
    )
    return reference_id, variant_sets


@pytest.fixture()
def gm():
    return GeneAminoAcidChangeToDNAVariants(
        reference=f"{DATA_DIR}/NC_000962.3.fasta",
        genbank=f"{DATA_DIR}/NC_000962.3.gb",
    )


def test_simple_gene(reference_seq):
    g = Gene(name="rpoB", reference=reference_seq, start=759807, end=763325)
    assert g.name == "rpoB"
    assert g.forward
    assert g.strand == "forward"
    assert (
        g.seq
        == "TTGGCAGATTCCCGCCAGAGCAAAACAGCCGCTAGTCCTAGTCCGAGTCGCCCGCAAAGTTCCTCGAATAACTCCGTACCCGGAGCGCCAAACCGGGTCTCCTTCGCTAAGCTGCGCGAACCACTTGAGGTTCCGGGACTCCTTGACGTCCAGACCGATTCGTTCGAGTGGCTGATCGGTTCGCCGCGCTGGCGCGAATCCGCCGCCGAGCGGGGTGATGTCAACCCAGTGGGTGGCCTGGAAGAGGTGCTCTACGAGCTGTCTCCGATCGAGGACTTCTCCGGGTCGATGTCGTTGTCGTTCTCTGACCCTCGTTTCGACGATGTCAAGGCACCCGTCGACGAGTGCAAAGACAAGGACATGACGTACGCGGCTCCACTGTTCGTCACCGCCGAGTTCATCAACAACAACACCGGTGAGATCAAGAGTCAGACGGTGTTCATGGGTGACTTCCCGATGATGACCGAGAAGGGCACGTTCATCATCAACGGGACCGAGCGTGTGGTGGTCAGCCAGCTGGTGCGGTCGCCCGGGGTGTACTTCGACGAGACCATTGACAAGTCCACCGACAAGACGCTGCACAGCGTCAAGGTGATCCCGAGCCGCGGCGCGTGGCTCGAGTTTGACGTCGACAAGCGCGACACCGTCGGCGTGCGCATCGACCGCAAACGCCGGCAACCGGTCACCGTGCTGCTCAAGGCGCTGGGCTGGACCAGCGAGCAGATTGTCGAGCGGTTCGGGTTCTCCGAGATCATGCGATCGACGCTGGAGAAGGACAACACCGTCGGCACCGACGAGGCGCTGTTGGACATCTACCGCAAGCTGCGTCCGGGCGAGCCCCCGACCAAAGAGTCAGCGCAGACGCTGTTGGAAAACTTGTTCTTCAAGGAGAAGCGCTACGACCTGGCCCGCGTCGGTCGCTATAAGGTCAACAAGAAGCTCGGGCTGCATGTCGGCGAGCCCATCACGTCGTCGACGCTGACCGAAGAAGACGTCGTGGCCACCATCGAATATCTGGTCCGCTTGCACGAGGGTCAGACCACGATGACCGTTCCGGGCGGCGTCGAGGTGCCGGTGGAAACCGACGACATCGACCACTTCGGCAACCGCCGCCTGCGTACGGTCGGCGAGCTGATCCAAAACCAGATCCGGGTCGGCATGTCGCGGATGGAGCGGGTGGTCCGGGAGCGGATGACCACCCAGGACGTGGAGGCGATCACACCGCAGACGTTGATCAACATCCGGCCGGTGGTCGCCGCGATCAAGGAGTTCTTCGGCACCAGCCAGCTGAGCCAATTCATGGACCAGAACAACCCGCTGTCGGGGTTGACCCACAAGCGCCGACTGTCGGCGCTGGGGCCCGGCGGTCTGTCACGTGAGCGTGCCGGGCTGGAGGTCCGCGACGTGCACCCGTCGCACTACGGCCGGATGTGCCCGATCGAAACCCCTGAGGGGCCCAACATCGGTCTGATCGGCTCGCTGTCGGTGTACGCGCGGGTCAACCCGTTCGGGTTCATCGAAACGCCGTACCGCAAGGTGGTCGACGGCGTGGTTAGCGACGAGATCGTGTACCTGACCGCCGACGAGGAGGACCGCCACGTGGTGGCACAGGCCAATTCGCCGATCGATGCGGACGGTCGCTTCGTCGAGCCGCGCGTGCTGGTCCGCCGCAAGGCGGGCGAGGTGGAGTACGTGCCCTCGTCTGAGGTGGACTACATGGACGTCTCGCCCCGCCAGATGGTGTCGGTGGCCACCGCGATGATTCCCTTCCTGGAGCACGACGACGCCAACCGTGCCCTCATGGGGGCAAACATGCAGCGCCAGGCGGTGCCGCTGGTCCGTAGCGAGGCCCCGCTGGTGGGCACCGGGATGGAGCTGCGCGCGGCGATCGACGCCGGCGACGTCGTCGTCGCCGAAGAAAGCGGCGTCATCGAGGAGGTGTCGGCCGACTACATCACTGTGATGCACGACAACGGCACCCGGCGTACCTACCGGATGCGCAAGTTTGCCCGGTCCAACCACGGCACTTGCGCCAACCAGTGCCCCATCGTGGACGCGGGCGACCGAGTCGAGGCCGGTCAGGTGATCGCCGACGGTCCCTGTACTGACGACGGCGAGATGGCGCTGGGCAAGAACCTGCTGGTGGCCATCATGCCGTGGGAGGGCCACAACTACGAGGACGCGATCATCCTGTCCAACCGCCTGGTCGAAGAGGACGTGCTCACCTCGATCCACATCGAGGAGCATGAGATCGATGCTCGCGACACCAAGCTGGGTGCGGAGGAGATCACCCGCGACATCCCGAACATCTCCGACGAGGTGCTCGCCGACCTGGATGAGCGGGGCATCGTGCGCATCGGTGCCGAGGTTCGCGACGGGGACATCCTGGTCGGCAAGGTCACCCCGAAGGGTGAGACCGAGCTGACGCCGGAGGAGCGGCTGCTGCGTGCCATCTTCGGTGAGAAGGCCCGCGAGGTGCGCGACACTTCGCTGAAGGTGCCGCACGGCGAATCCGGCAAGGTGATCGGCATTCGGGTGTTTTCCCGCGAGGACGAGGACGAGTTGCCGGCCGGTGTCAACGAGCTGGTGCGTGTGTATGTGGCTCAGAAACGCAAGATCTCCGACGGTGACAAGCTGGCCGGCCGGCACGGCAACAAGGGCGTGATCGGCAAGATCCTGCCGGTTGAGGACATGCCGTTCCTTGCCGACGGCACCCCGGTGGACATTATTTTGAACACCCACGGCGTGCCGCGACGGATGAACATCGGCCAGATTTTGGAGACCCACCTGGGTTGGTGTGCCCACAGCGGCTGGAAGGTCGACGCCGCCAAGGGGGTTCCGGACTGGGCCGCCAGGCTGCCCGACGAACTGCTCGAGGCGCAGCCGAACGCCATTGTGTCGACGCCGGTGTTCGACGGCGCCCAGGAGGCCGAGCTGCAGGGCCTGTTGTCGTGCACGCTGCCCAACCGCGACGGTGACGTGCTGGTCGACGCCGACGGCAAGGCCATGCTCTTCGACGGGCGCAGCGGCGAGCCGTTCCCGTACCCGGTCACGGTTGGCTACATGTACATCATGAAGCTGCACCACCTGGTGGACGACAAGATCCACGCCCGCTCCACCGGGCCGTACTCGATGATCACCCAGCAGCCGCTGGGCGGTAAGGCGCAGTTCGGTGGCCAGCGGTTCGGGGAGATGGAGTGCTGGGCCATGCAGGCCTACGGTGCTGCCTACACCCTGCAGGAGCTGTTGACCATCAAGTCCGATGACACCGTCGGCCGCGTCAAGGTGTACGAGGCGATCGTCAAGGGTGAGAACATCCCGGAGCCGGGCATCCCCGAGTCGTTCAAGGTGCTGCTCAAAGAACTGCAGTCGCTGTGCCTCAACGTCGAGGTGCTATCGAGTGACGGTGCGGCGATCGAACTGCGCGAAGGTGAGGACGAGGACCTGGAGCGGGCCGCGGCCAACCTGGGAATCAATCTGTCCCGCAACGAATCCGCAAGTGTCGAGGATCTTGCGTAA"
    )
    assert (
        g.prot
        == "LADSRQSKTAASPSPSRPQSSSNNSVPGAPNRVSFAKLREPLEVPGLLDVQTDSFEWLIGSPRWRESAAERGDVNPVGGLEEVLYELSPIEDFSGSMSLSFSDPRFDDVKAPVDECKDKDMTYAAPLFVTAEFINNNTGEIKSQTVFMGDFPMMTEKGTFIINGTERVVVSQLVRSPGVYFDETIDKSTDKTLHSVKVIPSRGAWLEFDVDKRDTVGVRIDRKRRQPVTVLLKALGWTSEQIVERFGFSEIMRSTLEKDNTVGTDEALLDIYRKLRPGEPPTKESAQTLLENLFFKEKRYDLARVGRYKVNKKLGLHVGEPITSSTLTEEDVVATIEYLVRLHEGQTTMTVPGGVEVPVETDDIDHFGNRRLRTVGELIQNQIRVGMSRMERVVRERMTTQDVEAITPQTLINIRPVVAAIKEFFGTSQLSQFMDQNNPLSGLTHKRRLSALGPGGLSRERAGLEVRDVHPSHYGRMCPIETPEGPNIGLIGSLSVYARVNPFGFIETPYRKVVDGVVSDEIVYLTADEEDRHVVAQANSPIDADGRFVEPRVLVRRKAGEVEYVPSSEVDYMDVSPRQMVSVATAMIPFLEHDDANRALMGANMQRQAVPLVRSEAPLVGTGMELRAAIDAGDVVVAEESGVIEEVSADYITVMHDNGTRRTYRMRKFARSNHGTCANQCPIVDAGDRVEAGQVIADGPCTDDGEMALGKNLLVAIMPWEGHNYEDAIILSNRLVEEDVLTSIHIEEHEIDARDTKLGAEEITRDIPNISDEVLADLDERGIVRIGAEVRDGDILVGKVTPKGETELTPEERLLRAIFGEKAREVRDTSLKVPHGESGKVIGIRVFSREDEDELPAGVNELVRVYVAQKRKISDGDKLAGRHGNKGVIGKILPVEDMPFLADGTPVDIILNTHGVPRRMNIGQILETHLGWCAHSGWKVDAAKGVPDWAARLPDELLEAQPNAIVSTPVFDGAQEAELQGLLSCTLPNRDGDVLVDADGKAMLFDGRSGEPFPYPVTVGYMYIMKLHHLVDDKIHARSTGPYSMITQQPLGGKAQFGGQRFGEMECWAMQAYGAAYTLQELLTIKSDDTVGRVKVYEAIVKGENIPEPGIPESFKVLLKELQSLCLNVEVLSSDGAAIELREGEDEDLERAAANLGINLSRNESASVEDLA"
    )


def test_reverse_gene(reference_seq):
    g = Gene(
        name="gidB",
        reference=reference_seq,
        start=4407528,
        end=4408202,
        forward=False,
    )
    assert g.name == "gidB"
    assert g.forward is False
    assert g.strand == "reverse"
    assert (
        g.seq
        == "ATGTCTCCGATCGAGCCCGCGGCGTCTGCGATCTTCGGACCGCGGCTTGGCCTTGCTCGGCGGTACGCCGAAGCGTTGGCGGGACCCGGTGTGGAGCGGGGGCTGGTGGGACCCCGCGAAGTCGGTAGGCTATGGGACCGGCATCTACTGAACTGCGCCGTGATCGGTGAGCTCCTCGAACGCGGTGACCGGGTCGTGGATATCGGTAGCGGAGCCGGGTTGCCGGGCGTGCCATTGGCGATAGCGCGGCCGGACCTCCAGGTAGTTCTCCTAGAACCGCTACTGCGCCGCACCGAGTTTCTTCGAGAGATGGTGACAGATCTGGGCGTGGCCGTTGAGATCGTGCGGGGGCGCGCCGAGGAGTCCTGGGTGCAGGACCAATTGGGCGGCAGCGACGCTGCGGTGTCACGGGCGGTGGCCGCGTTGGACAAGTTGACGAAATGGAGCATGCCGTTGATACGGCCGAACGGGCGAATGCTCGCCATCAAAGGCGAGCGGGCTCACGACGAAGTACGGGAGCACCGGCGTGTGATGATCGCATCGGGCGCGGTTGATGTCAGGGTGGTGACATGTGGCGCGAACTATTTGCGTCCGCCCGCGACCGTGGTGTTCGCACGACGTGGAAAGCAGATCGCCCGAGGGTCGGCACGGATGGCGAGTGGAGGGACGGCGTGA"
    )
    assert (
        g.prot
        == "MSPIEPAASAIFGPRLGLARRYAEALAGPGVERGLVGPREVGRLWDRHLLNCAVIGELLERGDRVVDIGSGAGLPGVPLAIARPDLQVVLLEPLLRRTEFLREMVTDLGVAVEIVRGRAEESWVQDQLGGSDAAVSRAVAALDKLTKWSMPLIRPNGRMLAIKGERAHDEVREHRRVMIASGAVDVRVVTCGANYLRPPATVVFARRGKQIARGSARMASGGTA"
    )


def test_reverse_gene2(reference_seq):
    g = Gene(
        name="katG",
        reference=reference_seq,
        start=2153889,
        end=2156111,
        forward=False,
    )
    assert g.name == "katG"
    assert g.forward is False
    assert g.strand == "reverse"
    assert (
        g.seq
        == "GTGCCCGAGCAACACCCACCCATTACAGAAACCACCACCGGAGCCGCTAGCAACGGCTGTCCCGTCGTGGGTCATATGAAATACCCCGTCGAGGGCGGCGGAAACCAGGACTGGTGGCCCAACCGGCTCAATCTGAAGGTACTGCACCAAAACCCGGCCGTCGCTGACCCGATGGGTGCGGCGTTCGACTATGCCGCGGAGGTCGCGACCATCGACGTTGACGCCCTGACGCGGGACATCGAGGAAGTGATGACCACCTCGCAGCCGTGGTGGCCCGCCGACTACGGCCACTACGGGCCGCTGTTTATCCGGATGGCGTGGCACGCTGCCGGCACCTACCGCATCCACGACGGCCGCGGCGGCGCCGGGGGCGGCATGCAGCGGTTCGCGCCGCTTAACAGCTGGCCCGACAACGCCAGCTTGGACAAGGCGCGCCGGCTGCTGTGGCCGGTCAAGAAGAAGTACGGCAAGAAGCTCTCATGGGCGGACCTGATTGTTTTCGCCGGCAACTGCGCGCTGGAATCGATGGGCTTCAAGACGTTCGGGTTCGGCTTCGGCCGGGTCGACCAGTGGGAGCCCGATGAGGTCTATTGGGGCAAGGAAGCCACCTGGCTCGGCGATGAGCGTTACAGCGGTAAGCGGGATCTGGAGAACCCGCTGGCCGCGGTGCAGATGGGGCTGATCTACGTGAACCCGGAGGGGCCGAACGGCAACCCGGACCCCATGGCCGCGGCGGTCGACATTCGCGAGACGTTTCGGCGCATGGCCATGAACGACGTCGAAACAGCGGCGCTGATCGTCGGCGGTCACACTTTCGGTAAGACCCATGGCGCCGGCCCGGCCGATCTGGTCGGCCCCGAACCCGAGGCTGCTCCGCTGGAGCAGATGGGCTTGGGCTGGAAGAGCTCGTATGGCACCGGAACCGGTAAGGACGCGATCACCAGCGGCATCGAGGTCGTATGGACGAACACCCCGACGAAATGGGACAACAGTTTCCTCGAGATCCTGTACGGCTACGAGTGGGAGCTGACGAAGAGCCCTGCTGGCGCTTGGCAATACACCGCCAAGGACGGCGCCGGTGCCGGCACCATCCCGGACCCGTTCGGCGGGCCAGGGCGCTCCCCGACGATGCTGGCCACTGACCTCTCGCTGCGGGTGGATCCGATCTATGAGCGGATCACGCGTCGCTGGCTGGAACACCCCGAGGAATTGGCCGACGAGTTCGCCAAGGCCTGGTACAAGCTGATCCACCGAGACATGGGTCCCGTTGCGAGATACCTTGGGCCGCTGGTCCCCAAGCAGACCCTGCTGTGGCAGGATCCGGTCCCTGCGGTCAGCCACGACCTCGTCGGCGAAGCCGAGATTGCCAGCCTTAAGAGCCAGATCCGGGCATCGGGATTGACTGTCTCACAGCTAGTTTCGACCGCATGGGCGGCGGCGTCGTCGTTCCGTGGTAGCGACAAGCGCGGCGGCGCCAACGGTGGTCGCATCCGCCTGCAGCCACAAGTCGGGTGGGAGGTCAACGACCCCGACGGGGATCTGCGCAAGGTCATTCGCACCCTGGAAGAGATCCAGGAGTCATTCAACTCCGCGGCGCCGGGGAACATCAAAGTGTCCTTCGCCGACCTCGTCGTGCTCGGTGGCTGTGCCGCCATAGAGAAAGCAGCAAAGGCGGCTGGCCACAACATCACGGTGCCCTTCACCCCGGGCCGCACGGATGCGTCGCAGGAACAAACCGACGTGGAATCCTTTGCCGTGCTGGAGCCCAAGGCAGATGGCTTCCGAAACTACCTCGGAAAGGGCAACCCGTTGCCGGCCGAGTACATGCTGCTCGACAAGGCGAACCTGCTTACGCTCAGTGCCCCTGAGATGACGGTGCTGGTAGGTGGCCTGCGCGTCCTCGGCGCAAACTACAAGCGCTTACCGCTGGGCGTGTTCACCGAGGCCTCCGAGTCACTGACCAACGACTTCTTCGTGAACCTGCTCGACATGGGTATCACCTGGGAGCCCTCGCCAGCAGATGACGGGACCTACCAGGGCAAGGATGGCAGTGGCAAGGTGAAGTGGACCGGCAGCCGCGTGGACCTGGTCTTCGGGTCCAACTCGGAGTTGCGGGCGCTTGTCGAGGTCTATGGCGCCGATGACGCGCAGCCGAAGTTCGTGCAGGACTTCGTCGCTGCCTGGGACAAGGTGATGAACCTCGACAGGTTCGACGTGCGCTGA"
    )
    assert (
        g.prot
        == "VPEQHPPITETTTGAASNGCPVVGHMKYPVEGGGNQDWWPNRLNLKVLHQNPAVADPMGAAFDYAAEVATIDVDALTRDIEEVMTTSQPWWPADYGHYGPLFIRMAWHAAGTYRIHDGRGGAGGGMQRFAPLNSWPDNASLDKARRLLWPVKKKYGKKLSWADLIVFAGNCALESMGFKTFGFGFGRVDQWEPDEVYWGKEATWLGDERYSGKRDLENPLAAVQMGLIYVNPEGPNGNPDPMAAAVDIRETFRRMAMNDVETAALIVGGHTFGKTHGAGPADLVGPEPEAAPLEQMGLGWKSSYGTGTGKDAITSGIEVVWTNTPTKWDNSFLEILYGYEWELTKSPAGAWQYTAKDGAGAGTIPDPFGGPGRSPTMLATDLSLRVDPIYERITRRWLEHPEELADEFAKAWYKLIHRDMGPVARYLGPLVPKQTLLWQDPVPAVSHDLVGEAEIASLKSQIRASGLTVSQLVSTAWAAASSFRGSDKRGGANGGRIRLQPQVGWEVNDPDGDLRKVIRTLEEIQESFNSAAPGNIKVSFADLVVLGGCAAIEKAAKAAGHNITVPFTPGRTDASQEQTDVESFAVLEPKADGFRNYLGKGNPLPAEYMLLDKANLLTLSAPEMTVLVGGLRVLGANYKRLPLGVFTEASESLTNDFFVNLLDMGITWEPSPADDGTYQGKDGSGKVKWTGSRVDLVFGSNSELRALVEVYGADDAQPKFVQDFVAAWDKVMNLDRFDVR"
    )


def test_get_codon(reference_seq):
    g = Gene(name="rpoB", reference=reference_seq, start=759807, end=763325)
    with pytest.raises(ValueError):
        g.get_codon(1173)
    assert g.get_codon(2) == "GCA"
    assert g.get_codon(3) == "GAT"
    assert g.get_reference_position(1) == 759807
    assert g.seq[0] == reference_seq[759806]
    assert g.get_reference_position(-1) == 759806


def test_get_codon_reverse(reference_seq):
    g = Gene(
        name="gidB",
        reference=reference_seq,
        start=4407528,
        end=4408202,
        forward=False,
    )
    with pytest.raises(ValueError):
        g.get_codon(225)
    assert g.get_codon(2) == "TCT"
    assert g.get_codon(3) == "CCG"
    assert g.get_reference_position(1) == 4408202
    assert g.get_reference_position(2) == 4408201
    assert g.get_reference_position(-1) == 4408203
    assert g.get_reference_position(-2) == 4408204


def test_gene_muts(gm):
    gm = GeneAminoAcidChangeToDNAVariants(
        reference=f"{DATA_DIR}/NC_000962.3.fasta",
        genbank=f"{DATA_DIR}/NC_000962.3.gb",
    )
    assert gm.get_alts("K") == ["AAA", "AAG"]
    # GAT -> ['GCA', 'GCT', 'GCC', 'GCG'], positions 759813,14,15
    assert sorted(gm.get_variant_names("rpoB", "D3A")) == sorted(
        ["GAT759813GCA", "GAT759813GCT", "GAT759813GCC", "GAT759813GCG"]
    )
    # GAT -> ['GCA', 'GCT', 'GCC', 'GCG'], positions 759813,14,15
    assert sorted(gm.get_variant_names("rpoB", "D3X")) == sorted(
        [
            "GAT759813GCA",
            "GAT759813GCT",
            "GAT759813GCC",
            "GAT759813GCG",
            "GAT759813TGT",
            "GAT759813TGC",
            "GAT759813GAA",
            "GAT759813GAG",
            "GAT759813GGA",
            "GAT759813GGT",
            "GAT759813GGC",
            "GAT759813GGG",
            "GAT759813TTT",
            "GAT759813TTC",
            "GAT759813ATA",
            "GAT759813ATT",
            "GAT759813ATC",
            "GAT759813CAT",
            "GAT759813CAC",
            "GAT759813AAA",
            "GAT759813AAG",
            "GAT759813ATG",
            "GAT759813TTA",
            "GAT759813TTG",
            "GAT759813CTA",
            "GAT759813CTT",
            "GAT759813CTC",
            "GAT759813CTG",
            "GAT759813AAT",
            "GAT759813AAC",
            "GAT759813CAA",
            "GAT759813CAG",
            "GAT759813CCA",
            "GAT759813CCT",
            "GAT759813CCC",
            "GAT759813CCG",
            "GAT759813AGT",
            "GAT759813AGC",
            "GAT759813TCA",
            "GAT759813TCT",
            "GAT759813TCC",
            "GAT759813TCG",
            "GAT759813AGA",
            "GAT759813AGG",
            "GAT759813CGA",
            "GAT759813CGT",
            "GAT759813CGC",
            "GAT759813CGG",
            "GAT759813ACA",
            "GAT759813ACT",
            "GAT759813ACC",
            "GAT759813ACG",
            "GAT759813TGG",
            "GAT759813GTA",
            "GAT759813GTT",
            "GAT759813GTC",
            "GAT759813GTG",
            "GAT759813TAT",
            "GAT759813TAC",
        ]
    )


def test_gene_muts2(gm):
    gm = GeneAminoAcidChangeToDNAVariants(
        reference=f"{DATA_DIR}/NC_000962.3.fasta",
        genbank=f"{DATA_DIR}/NC_000962.3.gb",
    )
    assert gm.get_alts("K") == ["AAA", "AAG"]
    # AGC -> ['CTT', 'CTC', 'CTA', 'CTG']
    #   # GAG -> ['GCA', 'GCT', 'GCC', 'GCG']
    # RC : CTC -> ['TGC',...] position2156103
    assert sorted(gm.get_variant_names("katG", "E3A")) == sorted(
        ["CTC2156103TGC", "CTC2156103AGC", "CTC2156103GGC", "CTC2156103CGC"]
    )


def test_make_variant_panel1(reference_id_and_variant_sets, gm):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta", kmer=31)
    gene = gm.get_gene("rpoB")
    for var in gm.get_variant_names("rpoB", "D3A"):
        ref, start, alt = split_var_name(var)
        v = Variant.create(
            variant_sets=variant_sets,
            reference=reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt],
        )
        panel = ag.create(v)
        for alt in panel.alts:
            seq = copy.copy(str(gene.seq))
            assert Seq(seq).translate()[2] == "D"
            seq = seq.replace(panel.refs[0][25:], alt[24:])
            assert seq != str(gene.seq)
            assert Seq(seq).translate()[2] == "A"


def test_make_variant_panel2(reference_id_and_variant_sets, gm):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta", kmer=31)
    gene = gm.get_gene("katG")
    for var in gm.get_variant_names("katG", "E3A"):
        ref, start, alt = split_var_name(var)
        v = Variant.create(
            variant_sets=variant_sets,
            reference=reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt],
        )
        panel = ag.create(v)
        for alt in panel.alts:
            seq = copy.copy(str(gene.seq.reverse_complement()))
            seq = seq.replace(
                panel.refs[0][:39], alt[: 39 + len(alt) - len(panel.refs[0])]
            )
            assert seq != str(gene.seq)
            # start with first codon in seq, so is multiple of 3, otherwise biopython warning
            assert (
                Seq(seq[-(len(gene.seq) - 3) :]).reverse_complement().translate()[2]
                == "A"
            )


def test_make_variant_panel3(reference_id_and_variant_sets, gm):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta")
    gene = gm.get_gene("katG")
    for var in gm.get_variant_names("katG", "S315L"):
        ref, start, alt = split_var_name(var)
        v = Variant.create(
            variant_sets=variant_sets,
            reference=reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt],
        )
        panel = ag.create(v)
        for alt in panel.alts:
            seq = copy.copy(str(gene.seq.reverse_complement()))
            seq = seq.replace(panel.refs[0], alt)
            assert seq != str(gene.seq)
            # start with first codon in seq, so is multiple of 3, otherwise biopython warning
            assert (
                Seq(seq[-(len(gene.seq) - 3) :]).reverse_complement().translate()[314]
                == "L"
            )


def test_make_variant_panel4(reference_id_and_variant_sets, gm):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta")
    gene = gm.get_gene("katG")
    for var in gm.get_variant_names("katG", "W90R"):
        ref, start, alt = split_var_name(var)
        v = Variant.create(
            variant_sets=variant_sets,
            reference=reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt],
        )
        panel = ag.create(v)
        for alt in panel.alts:
            seq = copy.copy(str(gene.seq.reverse_complement()))
            seq = seq.replace(panel.refs[0], alt)
            assert seq != str(gene.seq)
            # start with first codon in seq, so is multiple of 3, otherwise biopython warning
            assert (
                Seq(seq[-(len(gene.seq) - 3) :]).reverse_complement().translate()[89]
                == "R"
            )


def test_make_variant_panel5(reference_id_and_variant_sets, gm):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta")
    gene = gm.get_gene("gyrA")
    for var in gm.get_variant_names("gyrA", "D94X"):
        ref, start, alt = split_var_name(var)
        v = Variant.create(
            variant_sets=variant_sets,
            reference=reference_id,
            reference_bases=ref,
            start=start,
            alternate_bases=[alt],
        )
        panel = ag.create(v)
        for alt in panel.alts:
            seq = copy.copy(str(gene.seq))
            seq = seq.replace(panel.refs[0], alt)
            # start with first codon in seq, so is multiple of 3, otherwise biopython warning
            assert Seq(seq[-(len(gene.seq) - 3) :]).translate()[93] != "D"


def test_make_variant_panel6(reference_id_and_variant_sets, gm, reference_seq):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta", kmer=31)
    variants = list(gm.get_variant_names("pncA", "CAG28TAA", protein_coding_var=False))
    assert len(variants) == 1
    var = variants[0]
    ref, start, alt = split_var_name(var)
    assert ref == "CTG"
    assert start == 2289212
    assert alt == "TTA"
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference_id,
        reference_bases=ref,
        start=start,
        alternate_bases=[alt],
    )
    panel = ag.create(v)
    assert len(panel.alts) == 1
    alt = panel.alts[0]
    # the panel ref/alt seqs go past the end of the gene,
    # so can't comparie against gene sequence. Need to get
    # subseq from the reference seq
    panel_ref_start = reference_seq.find(panel.refs[0])
    assert panel_ref_start < start < panel_ref_start + len(panel.refs[0])
    seq = str(reference_seq[panel_ref_start : panel_ref_start + len(panel.refs[0])])
    assert seq == panel.refs[0]
    assert alt == seq[:30] + "TTA" + seq[33:]


def test_make_variant_panel7(reference_id_and_variant_sets, gm, reference_seq):
    reference_id, variant_sets = reference_id_and_variant_sets
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
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta", kmer=31)
    variants = list(gm.get_variant_names("eis", "G-10A", protein_coding_var=False))
    assert len(variants) == 1
    var = variants[0]
    ref, start, alt = split_var_name(var)
    assert ref == "C"
    assert start == 2715342
    assert alt == "T"
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference_id,
        reference_bases=ref,
        start=start,
        alternate_bases=[alt],
    )
    panel = ag.create(v)
    assert len(panel.alts) == 1
    alt = panel.alts[0]
    # the panel ref/alt seqs go past the end of the gene,
    # so can't comparie against gene sequence. Need to get
    # subseq from the reference seq
    panel_ref_start = reference_seq.find(panel.refs[0])
    assert panel_ref_start < start < panel_ref_start + len(panel.refs[0])
    seq = str(reference_seq[panel_ref_start : panel_ref_start + len(panel.refs[0])])
    assert seq == panel.refs[0]
    assert alt == seq[:30] + "T" + seq[31:]


def test_make_variant_panel8(reference_id_and_variant_sets, gm, reference_seq):
    reference_id, variant_sets = reference_id_and_variant_sets
    ag = AlleleGenerator(f"{DATA_DIR}/NC_000962.3.fasta", kmer=31)
    variants = list(gm.get_variant_names("eis", "TG-1T", protein_coding_var=False))
    assert len(variants) == 1
    var = variants[0]
    ref, start, alt = split_var_name(var)
    assert ref == "CA"
    assert start == 2715332
    assert alt == "A"
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference_id,
        reference_bases=ref,
        start=start,
        alternate_bases=[alt],
    )
    panel = ag.create(v)
    assert len(panel.alts) == 1
    alt = panel.alts[0]
    # the panel ref/alt seqs go past the end of the gene,
    # so can't comparie against gene sequence. Need to get
    # subseq from the reference seq
    panel_ref_start = reference_seq.find(panel.refs[0])
    assert panel_ref_start < start < panel_ref_start + len(panel.refs[0])
    seq = str(reference_seq[panel_ref_start : panel_ref_start + len(panel.refs[0])])
    assert seq == panel.refs[0]
    assert alt == seq[:30] + seq[31:]


def test_make_variant_panel_stop_codon(gm):
    variants = list(gm.get_variant_names("katG", "W90*", protein_coding_var=True))
    assert len(variants) == 3
    refs = sorted([split_var_name(v)[0] for v in variants])
    alts = sorted([split_var_name(v)[-1] for v in variants])
    var = variants[0]
    ref, start, alt = split_var_name(var)
    assert start == 2155842
    assert refs == ["CCA"] * 3
    assert alts == sorted(["TTA", "CTA", "TCA"])
