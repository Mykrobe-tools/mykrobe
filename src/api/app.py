import os
from flask import Flask
from flask import request

## Celery setup
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

def send_results(type, results, url):
    ## POST /samples/:id/result { type: "…", result: { … } }
    r = requests.post(url, json={'type': type, 'result' : results})

## Predictor
from api.analyses import run_predictor

@celery.task()
def predictor(file, sample_id):
    results=run_predictor(file, sample_id)
    url=os.path.join(ATLAS_API, "experiments", sample_id, "results")
    send_results("predictor", results, url)

## BIGSI
from api.analyses import BigsiTaskManager

BIGSI_DB_PATH=os.environ.get("BIGSI_DB_PATH", "dbpath") 
TB_REFERENCE_PATH=os.environ.get("TB_REFERENCE_PATH", "ref.fa") 
TB_GENBANK_PATH=os.environ.get("TB_GENBANK_PATH", "ref.gb") 
BIGSI_TM=BigsiTaskManager(BIGSI_DB_PATH,TB_REFERENCE_PATH,TB_GENBANK_PATH)

import hashlib
def _hash(w):
    w=w.encode('utf-8')
    h = hashlib.md5(w)
    return h.hexdigest()[:24]

@celery.task()
def bigsi(query_type, query):
    out={}
    results= {
        "sequence":BIGSI_TM.seq_query,
        "dna-variant":BIGSI_TM.dna_variant_query,
        "protein-variant":BIGSI_TM.protein_variant_query
    }[query_type](query)
    out["results"]=results
    out["query"]=query
    query_id=_hash(json.dumps(query))
    url=os.path.join(ATLAS_API, "queries", query_id, "results")    
    send_results("bigsi", out, url)

@app.route('/analyses', methods=["POST"])
def main():
    data=request.get_json()
    file = data.get('file', '')
    sample_id = data.get('sample_id', '')
    res=predictor.delay(file, sample_id)
    return json.dumps({"result":"success", "task_id":str(res)}), 200

@app.route('/search', methods=["POST"])
def search():
    data=request.get_json()
    t = data.get('type', '')
    query = data.get('query', '')
    res=bigsi.delay(t, query)
    return json.dumps({"result":"success", "task_id":str(res)}), 200    

## testing experiments requests /experiments/:sample_id/results
@app.route('/experiments/<sample_id>/results', methods=["POST"])
def results(sample_id):
    return request.data,200

@app.route('/queries/<query_id>/results', methods=["POST"])
def query_results(query_id):
    print(request.data,200)
    return request.data,200


