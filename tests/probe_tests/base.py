from mykrobe.utils import seq_to_kmers
def assert_no_overlapping_kmers(panel):
    for i in range(len(panel.refs)):
        kmer_intersection=set(seq_to_kmers(panel.refs[i], 31)).intersection(set(seq_to_kmers(panel.alts[i], 31)))
        if kmer_intersection:
            print(panel.refs)
            print(panel.alts)
            print(kmer_intersection)
        assert not kmer_intersection
