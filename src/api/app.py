import os
from flask import Flask
from flask import request

## Celery setup
from analyses import PredictorTaskManager
from analyses import BigsiTaskManager
from analyses import DistanceTaskManager

from celery import Celery
CELERY_BROKER_URL=os.environ.get("CELERY_BROKER_URL", 'redis://localhost:6379') 
DEFAULT_OUTDIR=os.environ.get("DEFAULT_OUTDIR", "./") 
ATLAS_API=os.environ.get("ATLAS_API", "https://api.atlas-prod.makeandship.com/") 
print(ATLAS_API)

def make_celery(app):
    celery = Celery(app.import_name, backend=app.config['CELERY_RESULT_BACKEND'],
                    broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task
    class ContextTask(TaskBase):
        abstract = True
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    celery.Task = ContextTask
    return celery

app = Flask(__name__)
app.config.update(
    CELERY_BROKER_URL=CELERY_BROKER_URL,
    CELERY_RESULT_BACKEND=CELERY_BROKER_URL
)
celery = make_celery(app)



import json
import requests
import logging
import http.client as http_client
http_client.HTTPConnection.debuglevel = 1

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)
requests_log = logging.getLogger("requests.packages.urllib3")
requests_log.setLevel(logging.DEBUG)
requests_log.propagate = True

def send_results(type, results, url, sub_type=None):
    ## POST /samples/:id/result { type: "…", result: { … } }
    d={'type': type, 'result' : results}
    if sub_type:
        d["subType"]=sub_type
    r = requests.post(url, json=d)

## Predictor

@celery.task()
def predictor_task(file, sample_id):
    results=PredictorTaskManager(DEFAULT_OUTDIR).run_predictor(file, sample_id)
    url=os.path.join(ATLAS_API, "experiments", sample_id, "results")
    send_results("predictor", results, url)

@app.route('/analyses', methods=["POST"])
def predictor():
    data=request.get_json()
    file = data.get('file', '')
    sample_id = data.get('sample_id', '')
    res=predictor_task.delay(file, sample_id)
    return json.dumps({"result":"success", "task_id":str(res)}), 200    

## BIGSI

BIGSI_DB_PATH=os.environ.get("BIGSI_DB_PATH", "dbpath") 
TB_REFERENCE_PATH=os.environ.get("TB_REFERENCE_PATH", "ref.fa") 
TB_GENBANK_PATH=os.environ.get("TB_GENBANK_PATH", "ref.gb") 

import hashlib
def _hash(w):
    w=w.encode('utf-8')
    h = hashlib.md5(w)
    return h.hexdigest()[:24]

@celery.task()
def bigsi(query_type, query):
    bigsi_tm=BigsiTaskManager(BIGSI_DB_PATH,TB_REFERENCE_PATH,TB_GENBANK_PATH)
    out={}
    results= {
        "sequence":bigsi_tm.seq_query,
        "dna-variant":bigsi_tm.dna_variant_query,
        "protein-variant":bigsi_tm.protein_variant_query
    }[query_type](query)
    out["results"]=results
    out["query"]=query
    query_id=_hash(json.dumps(query))
    url=os.path.join(ATLAS_API, "queries", query_id, "results")    
    send_results(query_type, out, url)

@app.route('/search', methods=["POST"])
def search():
    data=request.get_json()
    t = data.get('type', '')
    query = data.get('query', '')
    res=bigsi.delay(t, query)
    return json.dumps({"result":"success", "task_id":str(res)}), 200    

## nearest-neighbour neighbour distance

TREE_SAMPLES=[1,2,3]
@celery.task()
def distance_task(sample_id, distance_type):
    if distance_type == "all":
        results=DistanceTaskManager().distance(sample_id, sort=True)
    elif distance_type == "tree-distance":
        results=DistanceTaskManager().distance(sample_id, samples=TREE_SAMPLES, sort=True)        
    elif distance_type == "nearest-neighbour":
        results=DistanceTaskManager().distance(sample_id, limit=10, sort=True)
    else:
        raise TypeError("%s is not a valid query" % distance_type)
    url=os.path.join(ATLAS_API, "experiments", sample_id, "results")
    send_results("distance", results, url, sub_type=distance_type)

@app.route('/distance', methods=["POST"])
def distance():
    data=request.get_json()
    sample_id = data.get('sample_id', '')
    distance_type = data.get('distance_type', 'all') 
    assert distance_type in ["all", "tree-distance", "nearest-neighbour"]
    res=distance_task.delay(sample_id, distance_type)
    response= json.dumps({"result":"success", "task_id":str(res)}), 200     
    return response

## testing experiments requests /experiments/:sample_id/results
@app.route('/experiments/<sample_id>/results', methods=["POST"])
def results(sample_id):
    return request.data,200

@app.route('/queries/<query_id>/results', methods=["POST"])
def query_results(query_id):
    print(request.data,200)
    return request.data,200
