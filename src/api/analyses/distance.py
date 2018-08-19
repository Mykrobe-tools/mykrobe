import redis

REDIS=redis.StrictRedis(db=2)
SAMPLES_KEY="samples"
INTERMEDIATE_RESULT_EXPIRY=300

class DistanceTaskManager():

    def __init__(self, redis=REDIS, expiry=INTERMEDIATE_RESULT_EXPIRY):
        self.redis=redis   
        self.expiry=expiry             
        self.samples=self.__get_samples()

    def __get_samples(self):
        return {s.decode("utf-8") for s in self.redis.smembers(SAMPLES_KEY)}

    def __intermediate_key(self, s1, s2):
        return "_".join([str(s1),str(s2)])

    def _build_xor(self, primary_sample, samples):
        pipe=self.redis.pipeline()
        for secondary_sample in samples:
            k=self.__intermediate_key(primary_sample,secondary_sample)                    
            if secondary_sample != primary_sample:
                pipe.bitop("xor",k,primary_sample,secondary_sample)
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

    def distance(self, primary_sample):
        samples=self.__get_samples()
        self._build_xor(primary_sample, samples)
        return self._count_xor(primary_sample, samples)
