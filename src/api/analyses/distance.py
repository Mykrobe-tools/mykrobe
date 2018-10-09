import sys
import os
sys.path.append("../src/api/")

import json
import redis
import operator
from collections import OrderedDict

REDIS=redis.StrictRedis(db=2)
SAMPLES_KEY="samples"
INTERMEDIATE_RESULT_EXPIRY=300
def sort_and_filter_distance_dict(d, max_distance, limit):
    sorted_d = sorted(d.items(), key=operator.itemgetter(1))
    sorted_d = [x for x in sorted_d if x[1]<=max_distance]
    if limit:
        sorted_d=sorted_d[:limit]
    return OrderedDict(sorted_d)




class DistanceTaskManager():

    def __init__(self, redis=REDIS, expiry=INTERMEDIATE_RESULT_EXPIRY):
        self.redis=redis   
        self.expiry=expiry             
        self.samples=self.__get_samples()

    def __get_samples(self):
        return {s.decode("utf-8") for s in self.redis.smembers(SAMPLES_KEY)}

    def __intermediate_key(self, s1, s2):
        return "_".join([str(s1),"xor",str(s2)])

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
        samples=[s for s in samples if  s != primary_sample]
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

    def distance(self, primary_sample, max_distance=12, samples=None, limit=None, sort=True):
        if samples is None:
            samples=self.__get_samples()
        if limit is not None:
            sort = True
        self._build_xor(primary_sample, samples)
        distances = self._count_xor(primary_sample, samples)
        if sort:
            distances=sort_and_filter_distance_dict(distances, max_distance, limit)
        return distances

    ## Insert
    def _add_sample(self, sample_id):
        self.redis.sadd(SAMPLES_KEY, sample_id)

    def insert_from_json(json_path):
        with open(json_path, 'r') as inf:
            res=json.load(inf)
        return self.insert(res)

    def insert(self, res):
        for sample_id,data in res.items():
            self._add_sample(sample_id)
            genotypes = self._create_genotype_bitarray(data["genotypes"])
            passed_filter = self._create_filtered_bitarray(data["filtered"])
            filtered_genotypes=[x for i,x in enumerate(genotypes) if passed_filter[i]]
            self._insert_genotype_bitarray(filtered_genotypes, sample_id=sample_id)

    def _insert_genotype_bitarray(self, bitarray, sample_id):
        pipe = self.redis.pipeline()
        for i, j in enumerate(bitarray):
            pipe.setbit(self.__genotype_bitarray_key(sample_id), i, j)
        pipe.execute()    

    def _create_genotype_bitarray(self, sorted_calls):
        bitarray = [int(call > 1) for call in sorted_calls]
        return bitarray

    def _create_filtered_bitarray(self, sorted_calls):
        return sorted_calls    


