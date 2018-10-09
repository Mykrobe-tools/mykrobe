import logging
import os
import csv
import json

"""Adds variants to the database"""

logger = logging.getLogger(__name__)
import redis
r = redis.StrictRedis()


class BitarrayDistance():

    def __init__(self, redis=r):
        self.redis = redis

    def insert(self, json_path, sample_id):
        bitarray = self._create_genotype_bitarray(self.data["genotypes"])
        self._insert_genotype_bitarray(bitarray, sample=sample_id)
        bitarray = self._create_filtered_bitarray(self.data["filtered"])
        self._insert_filter_bitarray(bitarray, sample_id=sample_id)

    def _insert_genotype_bitarray(self, bitarray, sample_id):
        pipe = r.pipeline()
        for i, j in enumerate(bitarray):
            pipe.setbit("_".join([sample_id, "genotypes"]), i, j)
        pipe.execute()

    def _insert_genotype_bitarray(self, bitarray, sample_id):
        pipe = r.pipeline()
        for i, j in enumerate(bitarray):
            pipe.setbit("_".join([sample_id, "filters"]), i, j)
        pipe.execute()        

    def _create_genotype_bitarray(self, sorted_calls):
        bitarray = [int(call > 1) for call in sorted_calls]
        return bitarray

    def _create_filtered_bitarray(self, sorted_calls):
        return sorted_calls