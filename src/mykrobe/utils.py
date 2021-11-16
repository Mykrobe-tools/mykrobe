import os
import hashlib
import re
import json
from typing import TextIO

from Bio.Seq import Seq


def check_args(args):
    if args.db_name is None:
        args.db_name = os.environ.get("DB_NAME")
    if args.db_name is None:
        raise ValueError(
            "db_name needs to be set. Either run with --db_name :db_name or export DB_NAME=:db_name"
        )
    return args


def make_hash(s):
    return hashlib.sha256(s.encode("ascii", errors="ignore")).hexdigest()


def make_var_hash(ref, pos, alts):
    var = "".join([ref, str(pos), "/".join(alts)])
    return make_hash(var)


def split_var_name(name):
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/\*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


def unique(l):
    seen = set()
    return [x for x in l if x not in seen and not seen.add(x)]


def flatten(l):
    return [item for sublist in l for item in sublist]


def get_params(url):
    params = {}
    try:
        p_str = url.split("?")[1]
    except IndexError:
        return params
    p_str = p_str.split(" ")[0]
    p_str = p_str.split("&")
    for p in p_str:
        k, v = p.split("=")
        params[k] = v
    return params


def median(lst):
    if not lst:
        return 0
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if lstLen % 2:
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1]) / 2.0


def load_json(f):
    with open(f, "r") as infile:
        return json.load(infile)


def lazyprop(fn):
    attr_name = "_" + fn.__name__

    @property
    def _lazyprop(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazyprop


def seq_to_kmers(seq, kmer_size):
    for i in range(len(seq) - kmer_size + 1):
        yield seq[i : i + kmer_size]


def _x_mutation_fixed_var_name(var_name):
    """Takes mutation name from results base_json. If it is an "X" amino acid
    mutation, returns new name, where X is fixed with correct amino acid name.
    Otherwise returns None"""
    #                             prefix pos1     pos2
    #                              ||||| |||     |||||||
    # Example var_name with an X: "katG_S315X-GCT2155167GGT"
    #                                   |   | |||       |||
    #                                 aa1 aa2 codon1   codon2
    match = re.match(
        r"""(?P<prefix>.*_)(?P<aa1>[A-Z])(?P<pos1>[0-9]+)(?P<aa2>[A-Z])-(?P<codon1>[ACGT]{3})(?P<pos2>[0-9]+)(?P<codon2>[ACGT]{3})""",
        var_name,
    )
    if match is None or match.group("aa2") != "X":
        return None

    try:
        codon1_translate = str(Seq(match.group("codon1")).translate())
        codon1_rev_translate = str(
            Seq(match.group("codon1")).reverse_complement().translate()
        )
    except:
        return None

    # Check which strand the gene is on by checking if first codon or its
    # reverse complement match the amino acid. This is a safe way of
    # determining the strand because there is no codon where its translation
    # and the translation of its reverse complement results in the same
    # amino acid (this is checked in the tests of this repo)
    if codon1_translate == match.group("aa1"):
        try:
            new_aa = str(Seq(match.group("codon2")).translate())
        except:
            return None
    elif codon1_rev_translate == match.group("aa1"):
        try:
            new_aa = str(Seq(match.group("codon2")).reverse_complement().translate())
        except:
            return None
    else:
        return None

    return (
        match.group("prefix")
        + match.group("aa1")
        + match.group("pos1")
        + new_aa
        + "-"
        + match.group("codon1")
        + match.group("pos2")
        + match.group("codon2")
    )


def fix_amino_acid_X_variants_keys(dict_to_fix):
    """The way panels and variants work mean that the 'any amino acid' variants
    look like eg H445X, where X is any amino acid other than H. Users want to
    know the actual change. This function changes all keys in dict_to_fix that
    have those variants, changing the X to the actual amino acid."""
    keys_to_replace = {}
    for key in dict_to_fix:
        new_key = _x_mutation_fixed_var_name(key)
        if new_key is not None:
            if new_key in keys_to_replace or new_key in dict_to_fix:
                raise KeyError(
                    f"The 'X' amino acid in mutation {key}, resolves to {new_key}, "
                    f"which already exists in the panel"
                )
            keys_to_replace[key] = new_key

    for key, new_key in keys_to_replace.items():
        dict_to_fix[new_key] = dict_to_fix[key]
        del dict_to_fix[key]


def get_first_chrom_name(fp: TextIO) -> str:
    # make sure file pointer is at beginning of file
    fp.seek(0)
    header = fp.readline().rstrip()
    if not header.startswith(">"):
        raise ValueError(f"Expected fasta file, but it did not start with '>'")
    return header[1:].split()[0]
