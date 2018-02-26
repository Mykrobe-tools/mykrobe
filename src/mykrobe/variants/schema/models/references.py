from mongoengine import Document
from mongoengine import StringField
from mongoengine import IntField
from mongoengine import ListField
from mongoengine import BooleanField
from mongoengine import FloatField
from mongoengine import ReferenceField

from mykrobe.variants.schema.models.base import CreateAndSaveMixin


class Reference(Document, CreateAndSaveMixin):

    """
    A `Reference` is a canonical assembled contig, intended to act as a
    reference coordinate space for other genomic annotations. A single
    `Reference` might represent the human chromosome 1, for instance.

    `Reference`s are designed to be immutable.
    """
    name = StringField(unique=True, required=True)
    md5checksum = StringField(unique=True, required=True)
    reference_sets = ListField(ReferenceField("ReferenceSet"), required=True)

    length = IntField()
    source_accessions = ListField(StringField())
    sourceURI = StringField()
    # A sequence X is said to be derived from source sequence Y, if X and Y
    # are of the same length and the per-base sequence divergence at A/C/G/T bases
    # is sufficiently small. Two sequences derived from the same official
    # sequence share the same coordinates and annotations, and
    # can be replaced with the official sequence for certain use cases.
    is_derived = BooleanField()
    # The `sourceDivergence` is the fraction of non-indel bases that do not match the
    # reference this record was derived from.
    source_divergence = FloatField()
    # ID from http://www.ncbi.nlm.nih.gov/taxonomy (e.g. 9606->human)
    ncbi_taxon_id = IntField()

    @classmethod
    def create(cls, name, md5checksum, reference_sets,
               length=None, source_accessions=[],
               sourceURI=None, is_derived=None, source_divergence=None,
               ncbi_taxon_id=None):
        return cls(name=name, md5checksum=md5checksum,
                   reference_sets=reference_sets,
                   length=length, source_accessions=source_accessions,
                   sourceURI=sourceURI, is_derived=is_derived,
                   source_divergence=source_divergence,
                   ncbi_taxon_id=ncbi_taxon_id)


class ReferenceSet(Document, CreateAndSaveMixin):

    """
    A `ReferenceSet` is a set of `Reference`s which typically comprise a
    reference assembly, such as `GRCh38`. A `ReferenceSet` defines a common
    coordinate space for comparing reference-aligned experimental data.
    """

    """ The reference set name. """
    name = StringField(required=True)

    """ Optional free text description of this reference set. """
    description = StringField(default=None)

    # next information about the source of the sequences
    """ Public id of this reference set, such as `GRCh37`. """
    assembly_id = StringField()

    """ Specifies a FASTA format file/string. """
    sourceURI = StringField()

    """
    All known corresponding accession IDs in INSDC (GenBank/ENA/DDBJ) ideally
    with a version number, e.g. `NC_000001.11`.
    """
    sourceURI = ListField(StringField())

    """
    A reference set may be derived from a source if it contains
    additional sequences, or some of the sequences within it are derived
    (see the definition of `isDerived` in `Reference`).
    """
    is_derived = BooleanField()

    """
    ID from http://www.ncbi.nlm.nih.gov/taxonomy (e.g. 9606->human) indicating
    the species which this assembly is intended to model. Note that contained
    `Reference`s may specify a different `ncbiTaxonId`, as assemblies may
    contain reference sequences which do not belong to the modeled species, e.g.
    EBV in a human reference genome.
    """
    ncbi_taxon_id = IntField()

    @classmethod
    def create(cls, name, description=None, assembly_id=None,
               sourceURI=None, is_derived=None,
               ncbi_taxon_id=None):
        return cls(name=name, description=description,
                   assembly_id=assembly_id,
                   sourceURI=sourceURI, is_derived=is_derived,
                   ncbi_taxon_id=ncbi_taxon_id)

    @property
    def references(self):
        """
        References contained within this reference set
        """
        return Reference.objects(reference_sets=self)

    @property
    def md5checksum(self):
        """
        Order-independent MD5 checksum which identifies this `ReferenceSet`.

        To compute this checksum, make a list of `Reference.md5checksum` for all
        `Reference`s in this set. Then sort that list, and take the MD5 hash of
        all the strings concatenated together. Express the hash as a lower-case
        hexadecimal string.
        """
        raise NotImplementedError("")
