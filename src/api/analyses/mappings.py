import sys
import os

sys.path.append("../src/api/")
import redis

REDIS = redis.StrictRedis(db=3)
MAPPINGS_HKEY = "mappings"


def decode_dict(d):
    dd = {}
    for i, j in d.items():
        dd[i.decode("utf-8")] = j.decode("utf-8")
    return dd


class MappingsManager:
    def __init__(self, redis=REDIS, hkey=MAPPINGS_HKEY):
        self.hkey = hkey
        self.rhkey = "reverse" + self.hkey
        self.redis = redis

    def create_mapping(self, experiment_id, isolate_id):
        self.redis.hset(self.hkey, experiment_id, isolate_id)
        self.redis.hset(self.rhkey, isolate_id, experiment_id)

    def experiment_ids_to_isolate_ids(self):
        return decode_dict(self.redis.hgetall(self.hkey))

    def isolate_ids_to_experiment_ids(self):
        return decode_dict(self.redis.hgetall(self.rhkey))
