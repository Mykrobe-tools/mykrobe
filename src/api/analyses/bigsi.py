from bigsi.variants import BIGSIVariantSearch
from bigsi.variants import BIGSIAminoAcidMutationSearch
from bigsi.graph import BIGSI

DNA_QUERY_NAME="DNA"
PROT_QUERY_NAME="PROT"

class BigsiTaskManager():

	def __init__(self, bigsi_db_path, reference, genbank):
		self.bigsi=BIGSI(bigsi_db_path)
		self.bigsi_variant_search=BIGSIVariantSearch(self.bigsi, reference=reference)
		self.bigsi_protien_search=BIGSIAminoAcidMutationSearch(self.bigsi, reference=reference, genbank=genbank)


	def _filter_genotypes(self, results):
		o={}
		for var_name, _results in results.items():
			o[var_name]={}
			for s,v in _results.items():
				if sum([int(i) for i in v["genotype"].split("/") if i != "."])>=1:
					o[var_name][s]=v
		return o

	def seq_query(self, query):
		return self.bigsi.search(**query)

	def dna_variant_query(self, query):
		return self._filter_genotypes(self.bigsi_variant_search.search(**query))

	def protein_variant_query(self, query):
		return self._filter_genotypes(self.bigsi_protien_search.search(**query))


