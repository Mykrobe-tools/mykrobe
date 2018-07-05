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