import os
import hashlib
import re
import json

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
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/]+)", name, re.I).groups()
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
    #                              (prefix-) (--middle-)
    # Example var_name with an X: "katG_S315X-GCT2155167GGT"
    #                                       |           |||
    #                                       aa          codon
    match = re.match(
        r"""(?P<prefix>.*_[A-Z][0-9]+)(?P<aa>[A-Z])(?P<middle>-[ACGT]{3}[0-9]+)(?P<codon>[ACGT]{3})""",
        var_name,
    )
    if match is None or match.group("aa") != "X":
        return None
    try:
        amino_acid = str(Seq(match.group("codon")).translate())
    except:
        return None
    return (
        match.group("prefix")
        + amino_acid
        + match.group("middle")
        + match.group("codon")
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
            assert new_key not in keys_to_replace
            assert new_key not in dict_to_fix
            keys_to_replace[key] = new_key

    for key, new_key in keys_to_replace.items():
        dict_to_fix[new_key] = dict_to_fix[key]
        del dict_to_fix[key]
