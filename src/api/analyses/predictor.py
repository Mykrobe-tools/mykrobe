from api.helpers import load_json
from mykrobe.cmds.amr import run as predictor_cli
from mykrobe.parser import parser
import subprocess

def run_predictor(file, sample_id):
    cmd="mykrobe predict {sample_id} tb -1 {file} --output {sample_id}.json".format(sample_id=sample_id,file=file)
    outfile=os.path.join(DEFAULT_OUTDIR, "{sample_id}.json".format(sample_id=sample_id))
    out=subprocess.check_output(['mykrobe','predict', sample_id, "tb", "-1", file, "--output", outfile ])
    ## Load the output
    results=load_json(outfile)	