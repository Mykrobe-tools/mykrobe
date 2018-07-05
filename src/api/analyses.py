import os
from flask import Flask
from flask import request

from mykrobe.cmds.amr import run as predictor_cli
from mykrobe.parser import parser

## Celery setup
from celery import Celery

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

CELERY_BROKER_URL=os.environ.get("CELERY_BROKER_URL", 'redis://localhost:6379') 

app = Flask(__name__)
app.config.update(
    CELERY_BROKER_URL=CELERY_BROKER_URL,
    CELERY_RESULT_BACKEND=CELERY_BROKER_URL
)
celery = make_celery(app)

DEFAULT_OUTDIR=os.environ.get("DEFAULT_OUTDIR", "./") 
ATLAS_API=os.environ.get("ATLAS_API", "https://api.atlas-prod.makeandship.com/") 

from mykrobe.version import __version__
import json
import subprocess
import requests

def load_json(f):
    with open(f, 'r') as infile:
        return json.load(infile)

@celery.task()
def predictor(file, sample_id):
    cmd="mykrobe predict {sample_id} tb -1 {file} --output {sample_id}.json".format(sample_id=sample_id,file=file)
    outfile=os.path.join(DEFAULT_OUTDIR, "{sample_id}.json".format(sample_id=sample_id))
    out=subprocess.check_output(['mykrobe','predict', sample_id, "tb", "-1", file, "--output", outfile ])
    ## Load the output
    results=load_json(outfile)
    ## POST /samples/:id/result { type: "…", result: { … } }
    url=os.path.join(ATLAS_API, "experiments", sample_id, "results")
    r = requests.post(url, json={'type': "predictor", 'result' : results})

@app.route('/analyses', methods=["POST"])
def main():
    data=request.get_json()
    file = data.get('file', '')
    sample_id = data.get('sample_id', '')
    res=predictor.delay(file, sample_id)
    return json.dumps({"result":"success", "task_id":str(res)}), 200

## testing experiments requests /experiments/:sample_id/results
@app.route('/experiments/<sample_id>/results', methods=["POST"])
def results(sample_id):
    return request.data,200

