import datetime
import json

from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField
from mongoengine import FloatField
from mongoengine import DictField
from mongoengine import GenericReferenceField
from mongoengine import BooleanField
from mongoengine import queryset_manager

from mykrobe.variants.schema.models.base import CreateAndSaveMixin
from mykrobe.utils import split_var_name
from mykrobe.utils import make_var_hash

# Based on ga4gh Variant schema http://ga4gh.org/#/schemas feb 2016 with
# ocassional changes


class VariantSetMetadata(Document, CreateAndSaveMixin):
    key = StringField(unique_with="variant_set")
    value = StringField()
    type = StringField()
    # The number of values that can be included in a field described by this
    # metadata.
    number = IntField()
    description = StringField()
    info = DictField()
    variant_set = ReferenceField("VariantSet")

    @classmethod
    def create(
            cls,
            key,
            value,
            type,
            variant_set,
            number=None,
            description=None,
            info=None):
        return cls(key=key, value=value, type=type, description=description,
                   info=info, variant_set=variant_set)


class VariantSet(Document, CreateAndSaveMixin):

    """
    `Variant` and `CallSet` both belong to a `VariantSet`.
    `VariantSet` belongs to a `Dataset`.
    The variant set is equivalent to a VCF file.
    """
    name = StringField(required=True, unique=True)
    dataset = ReferenceField('Dataset')
    reference_set = ReferenceField('ReferenceSet', required=True)

    @classmethod
    def create(cls, name, reference_set, dataset=None):
        c = cls(name=name, reference_set=reference_set, dataset=dataset)
        return c.save()

    # Optional metadata associated with this variant set.
    # This array can be used to store information about the variant set, such as information found
    # in VCF header fields, that isn't already available in first class fields
    # such as "name".
    @property
    def metadata(self):
        return VariantSetMetadata.objects(variant_set=self)


class CallSet(Document, CreateAndSaveMixin):
    meta = {'allow_inheritance': True}

    """
         A `CallSet` is a collection of variant calls for a particular sample.
         It belongs to a `VariantSet`. This is simillar to one column in VCF.


    """
    name = StringField(required=True, default=None, unique=True)
    sample_id = StringField(required=True)
    created_at = DateTimeField(default=datetime.datetime.now)
    updated_at = DateTimeField(required=True, default=datetime.datetime.now)

    info = DictField()

    @classmethod
    def create(cls):
        raise NotImplementedError("Use VariantCallSet")


class VariantCallSet(CallSet):

    # Break from ga4gh schema -
    # When can a call set exist in multiple variant sets? If you have a set of
    # calls that you want to add to multiple variant sets. I think this demands
    # that a variant can exist in multiple variant sets, something not allowed
    # by ga4gh schema.

    variant_sets = ListField(ReferenceField('VariantSet'))

    @classmethod
    def create(cls, name, variant_sets, sample_id=None, info={}):
        c = cls(name=name, variant_sets=variant_sets, sample_id=sample_id,
                info=info)
        return c.save()


def convert_string_gt_to_list_int_gt(variant, genotype):
    try:
        return [int(i) for i in genotype.split('/')]
    except ValueError:
        return []


class Call(Document, CreateAndSaveMixin):
    meta = {'allow_inheritance': True,
            'indexes': [
                {
                    'fields': ['call_set']
                }
            ]
            }
    """
    A `Call` represents the determination of genotype with respect to a
    particular `Variant`.
    It may include associated information such as quality
    and phasing. For example, a call might assign a probability of 0.32 to
    the occurrence of a SNP named rs1234 in a call set with the name NA12345.
    """
    """
    The name of the call set this variant call belongs to.
    If this field is not present, the ordering of the call sets from a
    `SearchCallSetsRequest` over this `VariantSet` is guaranteed to match
    the ordering of the calls on this `Variant`.
    The number of results will also be the same.
    """
    call_set = ReferenceField('CallSet', required=True)
    """
    The genotype of this variant call.

    A 0 value represents the reference allele of the associated `Variant`. Any
    other value is a 1-based index into the alternate alleles of the associated
    `Variant`.

    If a variant had a referenceBases field of "T", an alternateBases
    value of ["A", "C"], and the genotype was [2, 1], that would mean the call
    represented the heterozygous value "CA" for this variant. If the genotype
    was instead [0, 1] the represented value would be "TA". Ordering of the
    genotype values is important if the phaseset field is present.
    """

    genotype = ListField(IntField())
    """
    The genotype likelihoods for this variant call. Each array entry
    represents how likely a specific genotype is for this call as
    log10(P(data | genotype)), analogous to the GL tag in the VCF spec. The
    value ordering is defined by the GL tag in the VCF spec (below)

    GL : genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible
    genotypes given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same
    ploidy is expected and the canonical order is used; without GT field, diploidy is assumed. If A is the allele in
    REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by:
    F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the
    ordering is: AA,AB,BB,AC,BC,CC, etc. For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Floats)

    """
    genotype_likelihoods = ListField(FloatField())

    info = DictField()

    @classmethod
    def create(cls):
        raise NotImplementedError("Use VariantCall")

    @property
    def call_set_name(self):
        return self.call_set.name

    @property
    def genotype_conf(self):
        # Returns the difference between the max and 2nd max liklihoods
        genotype_likelihoods = sorted(self.genotype_likelihoods, reverse=True)
        return genotype_likelihoods[0] - genotype_likelihoods[1]


class VariantCall(Call):

    meta = {'indexes': [
        {
            'fields': ['call_set']
        },
        {
            'fields': ['variant']
        }
    ]
    }
    variant = ReferenceField('Variant', required=True)  # Not in ga4gh
    # If this field is not null, this variant call's genotype ordering implies
    # the phase of the bases and is consistent with any other variant calls on
    # the same contig which have the same phaseset string.
    phaseset = GenericReferenceField(default=None)

    @classmethod
    def create(cls, variant, genotype, call_set=None, genotype_likelihoods=[],
               phaseset=None, info={}):
        if isinstance(genotype, str):
            genotype = convert_string_gt_to_list_int_gt(variant, genotype)
        if variant:
            cls._check_genotype_likelihood_length(
                genotype_likelihoods,
                variant)
        return cls(
            variant=variant,
            genotype=genotype,
            call_set=call_set,
            genotype_likelihoods=genotype_likelihoods,
            phaseset=phaseset,
            info=info)

    @classmethod
    def _check_genotype_likelihood_length(cls, genotype_likelihood, variant):
        if len(variant.alternate_bases) == 1:
            if not len(genotype_likelihood) == 3:
                raise ValueError(
                    "Biallelic sites should have 3 genotype likelihoods. AA,AB,BB")
        elif len(variant.alternate_bases) == 2:
            if not len(genotype_likelihood) == 6:
                raise ValueError(
                    "Biallelic sites should have 6 genotype likelihoods. AA,AB,BB,AC,BC,CC, etc")
        else:
            raise NotImplementedError(
                "Haven't implemented check for > triallelic sites")


class SequenceCall(Call):

    sequence = ReferenceField('Sequence', required=True)

    @classmethod
    def create(cls, sequence, call_set, genotype, genotype_likelihoods=[],
               info={}):
        if isinstance(genotype, str):
            genotype = convert_string_gt_to_list_int_gt(sequence, genotype)
        return cls(
            sequence=sequence,
            call_set=call_set,
            genotype=genotype,
            genotype_likelihoods=genotype_likelihoods,
            info=info)


def lazyprop(fn):
    attr_name = '_' + fn.__name__

    @property
    def _lazyprop(self):
        if not getattr(self, attr_name):
            setattr(self, attr_name, fn(self))
            self.save()
        return getattr(self, attr_name)
    return _lazyprop


def is_indel(reference_bases, alternate_bases):
    """ Return whether or not the variant is an INDEL """
    if len(reference_bases) > 1:
        return True
    for alt in alternate_bases:
        if alt is None:
            return True
        elif len(alt) != len(reference_bases):
            return True
    return False


def is_snp(reference_bases, alternate_bases):
    """ Return whether or not the variant is a SNP """
    if len(reference_bases) > 1:
        return False
    for alt in alternate_bases:
        if alt is None:
            return False
        if alt not in ['A', 'C', 'G', 'T', 'N', '*']:
            return False
    return True


def is_deletion(reference_bases, alternate_bases):
    """ Return whether or not the INDEL is a deletion """
    # if multiple alts, it is unclear if we have a transition
    if len(alternate_bases) > 1:
        return False

    if is_indel(reference_bases, alternate_bases):
        # just one alt allele
        alt_allele = alternate_bases[0]
        if alt_allele is None:
            return True
        if len(reference_bases) > len(alt_allele):
            return True
        else:
            return False
    else:
        return False


def is_insertion(reference_bases, alternate_bases):
    if len(alternate_bases) > 1:
        return False
    if is_indel(reference_bases, alternate_bases):
        # just one alt allele
        alt_allele = alternate_bases[0]
        if alt_allele is None:
            return False
        if len(alt_allele) > len(reference_bases):
            return True
        else:
            return False
    else:
        return False


def var_length(reference_bases, alternate_bases):
    return abs(len(reference_bases) - max([len(a) for a in alternate_bases]))


class Variant(Document, CreateAndSaveMixin):
    meta = {'indexes': [
        {
            'fields': ['start']
        },
        {
            'fields': ['var_hash']
        },
        {
            'fields': ['variant_sets']
        }
    ]
    }
    """A `Variant` represents a change in DNA sequence relative to some reference.
       For example, a variant could represent a SNP or an insertion.
      Variants belong to a `VariantSet`. This is simillar to a row in VCF.

      However, breaking from ga4gh we're allowing a variant belong to multiple
      VariantSets to allow a Variant to belong to multiple "Sets" of variants.
      """
    # Here, we've broken from ga4gh as they demand every variant belongs to a
    # single variant set. See CallSet.
    variant_sets = ListField(ReferenceField('VariantSet'), required=True)
    # The var_hash is a unique description on a variant. We use the hash of
    # "ref+pos+alt".
    var_hash = StringField(required=True)
    names = ListField(StringField())
    created_at = DateTimeField(required=True, default=datetime.datetime.now)
    updated_at = DateTimeField(required=True, default=datetime.datetime.now)
    start = IntField(required=True)  # (0-based)
    # The end position (exclusive), resulting in [start, end) closed-open
    # interval.
    end = IntField(required=False)
    reference_bases = StringField(required=True)
    alternate_bases = ListField(StringField(), required=True)
    info = DictField()

    # Each variant is defined against a single reference.
    # We can't have the same variant more than once i
    reference = ReferenceField(
        "Reference",
        required=True,
        unique_with="var_hash")

    length = IntField(required=False, default=None)
    is_snp = BooleanField(required=True)
    is_indel = BooleanField(required=True)
    is_deletion = BooleanField(required=True)
    is_insertion = BooleanField(required=True)

    @classmethod
    def create(cls, start, reference_bases,
               alternate_bases, variant_sets=None, reference=None, end=None,
               names=[], info={}):
        var_name = "".join(
            [reference_bases, str(start), "/".join(alternate_bases)])
        if var_name not in names:
            names.append(var_name)
        return cls(
            variant_sets=variant_sets,
            start=start,
            end=end,
            reference_bases=reference_bases,
            alternate_bases=alternate_bases,
            reference=reference,
            var_hash=make_var_hash(
                reference_bases,
                start,
                alternate_bases),
            names=names,
            info=info,
            length=var_length(reference_bases, alternate_bases),
            is_snp=is_snp(reference_bases, alternate_bases),
            is_indel=is_indel(reference_bases, alternate_bases),
            is_deletion=is_deletion(reference_bases, alternate_bases),
            is_insertion=is_insertion(reference_bases, alternate_bases)
        )

    @property
    def calls(self):
        # The variant calls for this particular variant. Each one represents the
        # determination of genotype with respect to this variant. `Call`s in this array
        # are implicitly associated with this `Variant`.
        return VariantCall.objects(variant=self)

    def count_calls(self):
        return VariantCall.objects(variant=self).count()

    @property
    def var_name(self):
        return "".join([self.reference_bases, str(
            self.start), "/".join(self.alternate_bases)])

    def add_to_variant_sets(self, variant_sets):
        self.update(add_to_set__variant_sets=variant_sets)

    def add_to_variant_set(self, variant_set):
        self.update(add_to_set__variant_sets=[variant_set])

    def __str__(self):
        return self.var_name

    def __repr__(self):
        return self.var_name

    def __eq__(self, other):
        return self.var_hash == other.var_hash

    def __gt__(self, other):
        return self.start > other.start

    def __lt__(self, other):
        return self.start < other.start

    @queryset_manager
    def snps(doc_cls, queryset):
        return queryset.filter(is_snp=True)

    @queryset_manager
    def indels(doc_cls, queryset):
        return queryset.filter(is_indel=True)

    @queryset_manager
    def deletions(doc_cls, queryset):
        return queryset.filter(is_deletion=True)

    @queryset_manager
    def insertions(doc_cls, queryset):
        return queryset.filter(is_insertion=True)

    @queryset_manager
    def ph_snps(doc_cls, queryset):
        return queryset.filter(
            is_indel=True,
            is_deletion=False,
            is_insertion=False)

    def overlapping(self, other):
        """Do these variants overlap in the reference"""
        return (
            other.start in self.ref_range) or (
            self.start in other.ref_range)

    @property
    def ref_range(self):
        return range(self.start, self.start + len(self.reference_bases))

    def split(self):
        if len(self.alternate_bases) == 1:
            return [self]
        variants = []
        for alt in self.alternate_bases:
            var = Variant.create(
                variant_sets=self.variant_sets,
                start=self.start,
                reference_bases=self.reference_bases,
                alternate_bases=[alt],
                reference=self.reference,
                end=self.end)
            variants.append(var)
        return variants

    def seen_in_samples(self):
        variant_call_sets = self.calls.distinct('call_set')
        return [variant_call_set.sample_id for variant_call_set in variant_call_sets]
