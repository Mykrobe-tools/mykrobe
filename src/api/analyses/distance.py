import redis
import operator
from collections import OrderedDict

REDIS=redis.StrictRedis(db=2)
SAMPLES_KEY="samples"
INTERMEDIATE_RESULT_EXPIRY=300
def sort_and_filter_distance_dict(d, limit):
    sorted_d = sorted(d.items(), key=operator.itemgetter(1))
    if limit:
        sorted_d=sorted_d[:limit]
    return OrderedDict(sorted_d)




class DistanceTaskManager():

    def __init__(self, redis=REDIS, expiry=INTERMEDIATE_RESULT_EXPIRY):
        self.inserter=BitarrayDistanceInserter(redis)
        self.redis=redis   
        self.expiry=expiry             
        self.samples=self.__get_samples()

    def __get_samples(self):
        return {s.decode("utf-8") for s in self.redis.smembers(SAMPLES_KEY)}

    def __intermediate_key(self, s1, s2):
        return "_".join([str(s1),str(s2)])

    def __genotype_bitarray_key(self,sample_id):
        return "_".join([sample_id, "genotypes"])        

    def _build_xor(self, primary_sample, samples):
        primary_sample_key=self.__genotype_bitarray_key(primary_sample)
        pipe=self.redis.pipeline()
        for secondary_sample in samples:
            secondary_sample_key=self.__genotype_bitarray_key(secondary_sample)
            k=self.__intermediate_key(primary_sample,secondary_sample)                    
            if secondary_sample != primary_sample:
                pipe.bitop("xor",k,primary_sample_key,secondary_sample_key)
                pipe.expire(k, self.expiry)
        pipe.execute()

    def _count_xor(self, primary_sample, samples):
        pipe=self.redis.pipeline()
        for secondary_sample in samples:
            k=self.__intermediate_key(primary_sample,secondary_sample)        
            if secondary_sample != primary_sample:
                pipe.bitcount(self.__intermediate_key(primary_sample,secondary_sample))
        res=pipe.execute()
        d={}
        for q,diff in zip(samples, res):
            d[q]=diff
        return d    

    def distance(self, primary_sample, samples=None, limit=None, sort=True):
        if samples is None:
            samples=self.__get_samples()
        if limit is not None:
            sort = True
        self._build_xor(primary_sample, samples)
        distances = self._count_xor(primary_sample, samples)
        if sort:
            distances=sort_and_filter_distance_dict(distances, limit)
        return distances

    ## Insert

    
    def insert(self, json_path, sample_id):
        genotypes = self._create_genotype_bitarray(self.data["genotypes"])
        filters = self._create_filtered_bitarray(self.data["filtered"])
        filtered_genotypes=[x for i,x in enumerate(genotypes) if not filters[i]]
        self._insert_genotype_bitarray(filtered_genotypes, sample=sample_id)

    def _insert_genotype_bitarray(self, bitarray, sample_id):
        pipe = r.pipeline()
        for i, j in enumerate(bitarray):
            pipe.setbit(self.__genotype_bitarray_key(sample_id), i, j)
        pipe.execute()    

    def _create_genotype_bitarray(self, sorted_calls):
        bitarray = [int(call > 1) for call in sorted_calls]
        return bitarray

    def _create_filtered_bitarray(self, sorted_calls):
        return sorted_calls    


