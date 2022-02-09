from __future__ import print_function

import os.path
import re
from collections import OrderedDict
from typing import Dict, Union, List, Tuple

import cyvcf2
from mongoengine import DoesNotExist
from mongoengine import NotUniqueError

from mykrobe import NON_METADATA_KEYS, SINGULAR_METADATA
from mykrobe.utils import make_var_hash
from mykrobe.variants.schema.models import ReferenceSet
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantCall
from mykrobe.variants.schema.models import VariantCallSet
from mykrobe.variants.schema.models import VariantSet
from mykrobe.variants.schema.models import VariantSetMetadata

GLOBAL_VARIANT_SET_NAME = "global_atlas"


class VCF(object):
    def __init__(
        self,
        f,
        reference_set_id,
        method="NotSpecified",
        force=False,
        append_to_global_variant_set=True,
    ):
        self.f = f
        self.reference_set = ReferenceSet.objects.get(id=reference_set_id)
        self.references = self._create_reference_lookup()
        self.method = method
        self.vcf_variant_set = None
        self.vcf_reader = None
        self.variants = []
        self.calls = []
        self.call_sets = {}
        self.force = force
        self.append_to_global_variant_set = append_to_global_variant_set

    def _create_reference_lookup(self):
        refs = {}
        for reference in self.reference_set.references:
            refs[reference.name] = reference
        return refs

    def add_to_database(self):
        self.vcf_reader = cyvcf2.VCF(self.f, gts012=True)
        self._create_new_variant_set()
        self._create_variant_set_meta_data()
        self._create_call_sets()
        self._create_variants_and_calls()
        self.vcf_reader.close()

    def _create_variants_and_calls(self):
        for record in self.vcf_reader:
            if not record.FILTER and self._is_record_valid(record):
                v = self._get_or_create_variant(record)
                for sample_idx in range(len(record.genotypes)):
                    gt_likelihoods = self._get_genotype_likelihoods(record, sample_idx)
                    c = VariantCall.create(
                        variant=v,
                        call_set=self.call_sets[self.vcf_reader.samples[sample_idx]],
                        genotype=str(Genotype(record.genotypes[sample_idx])),
                        genotype_likelihoods=gt_likelihoods,
                    )
                    self.calls.append(c)

        VariantCall.objects.insert(self.calls)

    def _get_or_create_variant(self, record):
        try:
            var_hash = make_var_hash(
                record.REF, record.POS, [str(a) for a in record.ALT]
            )
            v = Variant.objects.get(var_hash=var_hash)
            v.add_to_variant_set(self.vcf_variant_set)
        except DoesNotExist:
            try:
                reference = self.references[record.CHROM]
            except KeyError:
                raise KeyError(
                    "Reference %s cannot be found in reference set %s (%s). Please add it to the database."
                    % (record.CHROM, self.reference_set.id, self.reference_set.name)
                )
            v = Variant.create_and_save(
                variant_sets=self.variant_sets,
                start=record.POS,
                reference_bases=record.REF,
                alternate_bases=[str(a) for a in record.ALT],
                reference=reference,
                names=[record.ID],
            )
        return v

    def _remove_variant_set(self, variant_set_name):
        vs = VariantSet.objects.get(
            name=variant_set_name, reference_set=self.reference_set
        )
        for call_set in VariantCallSet.objects(variant_sets=vs):
            call_set.variant_sets.remove(vs)
            call_set.save()
            # Remove calls from callsets that only have this variantset
            if len(call_set.variant_sets) < 2:
                VariantCall.objects(call_set=call_set).delete()
                call_set.delete()
        # Remove variants that are ONLY from this variant set
        Variant.objects(variant_sets=vs, variant_sets__size=2).delete()
        VariantSetMetadata.objects(variant_set=vs).delete()
        vs.delete()

    def _create_new_variant_set(self):
        variant_set_name = os.path.basename(self.f)
        if VariantSet.objects(name=variant_set_name, reference_set=self.reference_set):
            if not self.force:
                raise NotUniqueError(
                    "VariantSet %s already exists. Rerun with -f to recreate."
                    % variant_set_name
                )
            else:
                self._remove_variant_set(variant_set_name)
        self.vcf_variant_set = VariantSet.create_and_save(
            name=variant_set_name, reference_set=self.reference_set
        )

    def _create_call_sets(self):
        for sample in self.vcf_reader.samples:
            try:
                cs = VariantCallSet.create_and_save(
                    name="_".join([sample, self.method]),
                    variant_sets=self.variant_sets,
                    sample_id=sample,
                    info={"variant_caller": self.method},
                )
            except NotUniqueError:
                raise ValueError(
                    "There is already a call set for sample %s with method %s "
                    % (sample, self.method)
                )
            else:
                self.call_sets[sample] = cs

    @property
    def variant_sets(self):
        if self.append_to_global_variant_set:
            return [self.vcf_variant_set, self.global_variant_set]
        else:
            return [self.vcf_variant_set]

    @property
    def global_variant_set(self):
        try:
            vs = VariantSet.objects.get(name=GLOBAL_VARIANT_SET_NAME)
        except:
            vs = VariantSet.create_and_save(
                name=GLOBAL_VARIANT_SET_NAME, reference_set=self.reference_set
            )
        return vs

    def _create_variant_set_meta_data(self):
        for variant_set in self.variant_sets:
            for k, v in self._metadata_from_header(self.vcf_reader.raw_header).items():
                if not VariantSetMetadata.objects(key=k, variant_set=variant_set):
                    vsm = VariantSetMetadata.create_and_save(
                        key=k,
                        value="metadata",
                        type="metadata",
                        variant_set=variant_set,
                    )
            for hrec in self.vcf_reader.header_iter():
                try:
                    hrec_id = hrec["ID"]
                    hdr_type = hrec["HeaderType"].upper()
                except KeyError:  # cyvcf2's HREC doesn't have a `get` method
                    continue

                if not VariantSetMetadata.objects(key=hrec_id, variant_set=variant_set):
                    try:
                        desc = hrec["Description"]
                    except KeyError:
                        desc = ""

                    params = dict(
                        key=hrec_id, variant_set=variant_set, description=desc
                    )
                    if hdr_type in ("INFO", "FORMAT"):
                        params["value"] = hdr_type.lower() + "s"
                        for field in ("type", "number"):
                            try:
                                params[field] = hrec[field.capitalize()]
                            except KeyError:
                                raise KeyError(
                                    "{} ID {} does not have field '{}'".format(
                                        hdr_type, hrec_id, field.capitalize()
                                    )
                                )
                    elif hdr_type == "FILTER":
                        params["value"] = "filters"
                        params["type"] = "filters"
                    else:
                        continue

                    vsm = VariantSetMetadata.create_and_save(**params)

    @staticmethod
    def _read_meta_hash(meta_string: str) -> Tuple[str, Dict[str, str]]:
        """Taken from pyvcf
        https://github.com/jamescasbon/PyVCF/blob/476169cd457ba0caa6b998b301a4d91e975251d9/vcf/parser.py#L184-L220
        """
        items = meta_string.split("=", 1)
        # Removing initial hash marks
        key = items[0].lstrip("#")
        # N.B., items can have quoted values, so cannot just split on comma
        val = OrderedDict()
        state = 0
        k = ""
        v = ""
        for c in items[1].strip("[<>]"):

            if state == 0:  # reading item key
                if c == "=":
                    state = 1  # end of key, start reading value
                else:
                    k += c  # extend key
            elif state == 1:  # reading item value
                if v == "" and c == '"':
                    v += c  # include quote mark in value
                    state = 2  # start reading quoted value
                elif c == ",":
                    val[k] = v  # store parsed item
                    state = 0  # read next key
                    k = ""
                    v = ""
                else:
                    v += c
            elif state == 2:  # reading quoted item value
                if c == '"':
                    v += c  # include quote mark in value
                    state = 1  # end quoting
                else:
                    v += c
        if k != "":
            val[k] = v
        return key, val

    def _read_meta(self, meta_string: str) -> Tuple[str, Union[str, Dict[str, str]]]:
        """Taken from pyvcf
        https://github.com/jamescasbon/PyVCF/blob/476169cd457ba0caa6b998b301a4d91e975251d9/vcf/parser.py#L222-L231
        """
        if re.match("##.+=<", meta_string):
            return self._read_meta_hash(meta_string)
        meta_pattern = re.compile(r"""##(?P<key>.+?)=(?P<val>.+)""")
        match = meta_pattern.match(meta_string)
        if not match:
            # Spec only allows key=value, but we try to be liberal and
            # interpret anything else as key=none (and all values are parsed
            # as strings).
            return meta_string.lstrip("#"), "none"
        return match.group("key"), match.group("val")

    def _metadata_from_header(
        self, raw_header: str
    ) -> Dict[str, Union[str, List[str]]]:
        metadata = dict()
        for line in raw_header.splitlines():
            key, val = self._read_meta(line)
            if key in NON_METADATA_KEYS:
                continue
            if key in SINGULAR_METADATA:
                metadata[key] = val
            else:
                if key not in metadata:
                    metadata[key] = []
                metadata[key].append(val)

        return metadata

    @staticmethod
    def _is_record_valid(record: cyvcf2.Variant) -> bool:
        valid = True
        for i, gt in enumerate(map(Genotype, record.genotypes)):
            if not gt.is_hom_alt():
                valid = False
                break

            gt_conf = record.format("GT_CONF")[i][0]
            if gt_conf <= 1:  # todo: magic number.
                valid = False
                break

        return valid

    @staticmethod
    def _get_genotype_likelihoods(
        record: cyvcf2.Variant, sample_idx: int
    ) -> List[float]:
        try:
            genotype_likelihoods = record.format("GL")[sample_idx]
        except KeyError:
            gt_conf = record.format("GT_CONF")[sample_idx][0]
            genotype_likelihoods = [0, 0, 0]
            genotype_likelihoods[sum(record.genotypes[sample_idx][:-1])] = gt_conf

        return genotype_likelihoods


class Genotype(object):
    """This small class makes it easier to go from a cyvcf2 genotype array to a string"""

    __slots__ = ("alleles", "phased")

    def __init__(self, genotype):
        """genotype is a numpy array whose first element is an int, and last element
        is a bool. There are a variable number of ints. The actual type is
        List[int, ..., bool], but this doesn't seem to be supported
        https://stackoverflow.com/q/71042017/5299417
        """
        self.alleles = genotype[:-1]
        self.phased = genotype[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        # return sep.join(["." if a == -1 else str(a) for a in self.alleles])
        return sep.join(map(self._allele_to_str, self.alleles))
        # return sep.join("0123."[a] for a in self.alleles)

    @staticmethod
    def _allele_to_str(a: int) -> str:
        return "." if a == -1 else str(a)

    def is_het(self) -> bool:
        return len(set(self.alleles)) > 1

    def is_null(self) -> bool:
        return sum(self.alleles) < 0

    def is_hom_alt(self) -> bool:
        """Is genotype homozygous alternate?"""
        return not self.is_het() and self.alleles[0] > 0

    __repr__ = __str__
