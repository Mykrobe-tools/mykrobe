# A simple API to trigger mykrobe-atlas-cli analyses

working dir /src/api
```
celery -A app.celery worker -E
ATLAS_API="localhost:8080" DEFAULT_OUTDIR="/atlas/predictor-results/" CELERY_BROKER_URL='redis://localhost:6379' FLASK_DEBUG=1 FLASK_APP=app.py flask run --port 8080
```

```
curl -H "Content-Type: application/json" -X POST -d '{"file":"path/to/file", "experiment_id": "experiment_id"}' localhost:8080/analyses
```

On jessie:

export PATH=/ssd0/software/mccortex/bin/:$PATH
curl -H "Content-Type: application/json" -X POST -d '{"file":"/atlas/test-data/MDR.fastq.gz", "experiment_id": "MDR_test"}' localhost:8080/analyses


## Search query
```
export TB_REFERENCE_PATH="../mykrobe/data/NC_000962.3.fasta" 
export TB_GENBANK_PATH="../mykrobe/data/NC_000962.3.gb" 
export BIGSI_DB_PATH="~/git/BIGSI/db-bigsi/"

export TB_REFERENCE_PATH="/home/admin/git/mykrobe-atlas-cli/src/mykrobe/data/NC_000962.3.fasta" 
export TB_GENBANK_PATH="/home/admin/git/mykrobe-atlas-cli/src/mykrobe/data/NC_000962.3.gb" 
export BIGSI_DB_PATH="/ssd1/sra/"

ATLAS_API="http://localhost:8080" celery -A app.celery worker -E


ATLAS_API="http://localhost:8080" DEFAULT_OUTDIR="/atlas/predictor-results/" CELERY_BROKER_URL='redis://localhost:6379' FLASK_DEBUG=1 FLASK_APP=app.py flask run --port 8080
```

```
curl -H "Content-Type: application/json" -X POST -d '{"type":"sequence","query":{"seq":"CGGTCAGTCCGTTTGTTCTTGTGGCGAGTGTTGCCGTTTTCTTG", "threshold":0.9, "user_id": "1234567", "result_id": "2345678" } }' localhost:8080/search
```

```
curl -H "Content-Type: application/json" -X POST -d '{"type":"dna-variant", "query": {"ref":"T","alt":"A","pos":99} }' localhost:8080/search
```
```
curl -H "Content-Type: application/json" -X POST -d '{"type":"protein-variant", "query": {"ref":"S","alt":"L","pos":450,"gene":"rpoB"} }' localhost:8080/search
```

## Distance
```
curl -H "Content-Type: application/json" -X POST -d '{"experiment_id": "experiment_id"}' localhost:8080/distance

curl -H "Content-Type: application/json" -X POST -d '{"experiment_id": "experiment_id", "distance_type":"tree-distance"}' localhost:8080/distance

curl -H "Content-Type: application/json" -X POST -d '{"experiment_id": "experiment_id", "distance_type":"nearest-neighbour"}' localhost:8080/distance

```

## Tree
```
curl localhost:8080/tree/latest
```
