from mongoengine import Document
from mongoengine import DictField

class MykrobePredictorPhylogeneticsResult(Document):

    phylogenetics = DictField()

    @classmethod
    def create(cls, phylogenetics):
        return cls(phylogenetics = phylogenetics)  

    def to_dict(self):
    	return self.to_mongo().to_dict()          
