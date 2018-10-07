import os
from flask import Flask
from flask import request
from Bio import Phylo
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
## Celery setup
from analyses import PredictorTaskManager
from analyses import BigsiTaskManager
from analyses import DistanceTaskManager
from analyses import MappingsManager

from celery import Celery
CELERY_BROKER_URL=os.environ.get("CELERY_BROKER_URL", 'redis://localhost:6379') 
DEFAULT_OUTDIR=os.environ.get("DEFAULT_OUTDIR", "./") 
ATLAS_API=os.environ.get("ATLAS_API", "https://api.atlas-prod.makeandship.com/") 
MAPPER=MappingsManager()
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

def send_results(type, results, url, sub_type=None, request_type="POST"):
    ## POST /isolates/:id/result { type: "…", result: { … } }
    d={'type': type, 'result' : results}
    if sub_type:
        d["subType"]=sub_type
    if request_type == "PUT":
        r=requests.put(url, json=d)
    else:
        r = requests.post(url, json=d)

## Predictor

@celery.task()
def predictor_task(file, experiment_id):
    results=PredictorTaskManager(DEFAULT_OUTDIR).run_predictor(file, experiment_id)
    url=os.path.join(ATLAS_API, "experiments", experiment_id, "results")
    send_results("predictor", results, url)

@celery.task()
def genotype_task(file, experiment_id):
    results=PredictorTaskManager(DEFAULT_OUTDIR).run_genotype(file, experiment_id)
    url=os.path.join(ATLAS_API, "experiments", experiment_id, "results")
    # send_results("genotype", results, url)
    ## Insert distance results 
    DistanceTaskManager().insert(results)
    # ## Trigger distance tasks
    res=distance_task.delay(experiment_id, "tree-distance")
    res=distance_task.delay(experiment_id, "nearest-neighbour")

@app.route('/analyses', methods=["POST"])
def predictor():
    data=request.get_json()
    file = data.get('file', '')
    experiment_id = data.get('experiment_id', '')
    res=predictor_task.delay(file, experiment_id)
    res=genotype_task.delay(file, experiment_id)
    MAPPER.create_mapping(experiment_id, experiment_id)
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
    user_id=query["user_id"]
    result_id=query["result_id"]
    url=os.path.join(ATLAS_API, "users", user_id, "results", result_id)    
    send_results(query_type, out, url, request_type="PUT")

@app.route('/search', methods=["POST"])
def search():
    data=request.get_json()
    t = data.get('type', '')
    query = data.get('query', '')
    res=bigsi.delay(t, query)
    return json.dumps({"result":"success", "task_id":str(res)}), 200    

## nearest-neighbour neighbour distance

## Tree
def load_tree(version):
    if version == "latest":
        float_versions=[float(x) for x in sorted(TREE_PATH.keys())]
        version_float_max=max(float_versions)
        float_versions_index=[i for i,x in enumerate(float_versions) if x==version_float_max]
        version=list(sorted(TREE_PATH.keys()))[float_versions_index[0]]
        tree_path=TREE_PATH[version]
    with open(tree_path, 'r') as infile:
        data=infile.read().replace('\n', '')    
    return data

@celery.task()
def tree_task(version):
    assert version == "latest"
    data=load_tree(version)
    url=os.path.join(ATLAS_API, "trees")
    results={"tree":data, "version":version}
    send_results("tree", results, url)
    return results

@app.route('/tree/<version>', methods=["GET"])
def tree(version):
    results=tree_task(version)
    response= json.dumps({"result":results, "type": "tree"}), 200 
    return response

TREE_PATH={"1.0":"data/tb_newick.txt"}
def get_tree_isolates():
    newick=load_tree("latest")
    tree=Phylo.read(StringIO(newick), 'newick')
    tree_isolates=[c.name for c in tree.get_terminals()]
    return tree_isolates


TREE_ISOLATES=get_tree_isolates()
DEFAULT_MAX_NN_DISTANCE=12
DEFAULT_MAX_NN_EXPERIMENTS=100
@celery.task()
def distance_task(experiment_id, distance_type, max_distance=None, limit=None):
    if max_distance is None:
        max_distance=DEFAULT_MAX_NN_DISTANCE
    if limit is None:
        limit=DEFAULT_MAX_NN_EXPERIMENTS
    if distance_type == "all":
        results=DistanceTaskManager().distance(experiment_id, sort=True)
    elif distance_type == "tree-distance":
        results=DistanceTaskManager().distance(experiment_id, isolates=TREE_ISOLATES, sort=True)        
    elif distance_type == "nearest-neighbour":
        results=DistanceTaskManager().distance(experiment_id, max_distance=max_distance, limit=limit, sort=True)
    else:
        raise TypeError("%s is not a valid query" % distance_type)
    url=os.path.join(ATLAS_API, "experiments", experiment_id, "results")
    send_results("distance", results, url, sub_type=distance_type)

@app.route('/distance', methods=["POST"])
def distance():
    data=request.get_json()
    experiment_id = data.get('experiment_id', '')
    distance_type = data.get('distance_type', 'all') 
    max_distance = data.get('max_distance') 
    limit = data.get('limit') 
    assert distance_type in ["all", "tree-distance", "nearest-neighbour"]
    res=distance_task.delay(experiment_id, distance_type, max_distance=max_distance, limit=limit)
    response= json.dumps({"result":"success", "task_id":str(res)}), 200     
    return response

## Mappings from experiment_id to isolate_id
@app.route('/mappings', methods=["GET"])
def mappings():
    mappings=MAPPER.experiment_ids_to_isolate_ids()
    return json.dumps(mappings)

@app.route('/reverse_mappings', methods=["GET"])
def mappings():
    mappings=MAPPER.isolate_ids_to_experiment_ids()
    return json.dumps(mappings)


## testing experiments requests /experiments/:experiment_id/results
@app.route('/experiments/<experiment_id>/results', methods=["POST"])
def results(experiment_id):
    return request.data,200

@app.route('/queries/<query_id>/results', methods=["POST"])
def query_results(query_id):
    return request.data,200

@app.route('/trees', methods=["POST"])
def tree_results():
    print(request.data)
    return request.data,200    
