A simple API to trigger mykrobe-atlas-cli analyses

working dir /src/api
```
celery -A analyses.celery worker -E
ATLAS_API="localhost:8080" DEFAULT_OUTDIR="/atlas/predictor-results/" CELERY_BROKER_URL='redis://localhost:6379' FLASK_DEBUG=1 FLASK_APP=analyses.py flask run --port 8080
```

```
curl -H "Content-Type: application/json" -X POST -d '{"file":"path/to/file", sample_id: "sample_id"}' localhost:8080/analyses
```

On jessie:

export PATH=/ssd0/software/mccortex/bin/:$PATH
curl -H "Content-Type: application/json" -X POST -d '{"file":"/atlas/test-data/MDR.fastq.gz", "sample_id": "MDR_test"}' localhost:8080/analyses


## Search query
```
export TB_REFERENCE_PATH="/Users/phelimb/git/mykrobe-atlas-cli/src/mykrobe/data/NC_000962.3.fasta" 
export TB_GENBANK_PATH="/Users/phelimb/git/mykrobe-atlas-cli/src/mykrobe/data/NC_000962.3.gb" 
export BIGSI_DB_PATH="/Users/phelimb/git/BIGSI/db-bigsi"

celery -A analyses.celery worker -E


ATLAS_API="localhost:8080" DEFAULT_OUTDIR="/atlas/predictor-results/" CELERY_BROKER_URL='redis://localhost:6379' FLASK_DEBUG=1 FLASK_APP=analyses.py flask run --port 8080
```

```
 curl -H "Content-Type: application/json" -X POST -d '{"type":"sequence", "query": {"seq":"ATCGATCG", "threshold":0.9} }' localhost:8080/search
```
```
curl -H "Content-Type: application/json" -X POST -d '{"type":"DNA-variant", query: {ref:"A",alt:"T","pos":99 "threshold":.9} }' localhost:8080/search
```
```
curl -H "Content-Type: application/json" -X POST -d '{"type":"protein-variant", query: {ref:"L",alt:"X","pos":450,gene:"rpoB","threshold":.9} }' localhost:8080/search
```