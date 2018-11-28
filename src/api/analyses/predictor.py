from helpers import load_json
import subprocess
import os


class PredictorTaskManager:
    def __init__(self, outdir):
        self.outdir = outdir

    def predictor_filepath(self, sample_id):
        return os.path.join(
            self.outdir, "{sample_id}_predictor.json".format(sample_id=sample_id)
        )

    def genotype_filepath(self, sample_id):
        return os.path.join(
            self.outdir, "{sample_id}_genotype.json".format(sample_id=sample_id)
        )

    def run_predictor(self, file, sample_id):
        outfile = self.predictor_filepath(sample_id)
        cmd = "mykrobe predict --keep_tmp {sample_id} tb -1 {file} --output {sample_id}.json".format(
            sample_id=sample_id, file=file
        )
        out = subprocess.check_output(
            [
                "mykrobe",
                "predict",
                sample_id,
                "tb",
                "--panel",
                "atlas",
                "-1",
                file,
                "--output",
                outfile,
            ]
        )
        ## Load the output
        results = load_json(outfile)
        return results

    def run_genotype(self, file, sample_id):
        outfile = self.genotype_filepath(sample_id)
        cmd = "mykrobe genotype --keep_tmp {sample_id} data/tb-k21-probe-set-feb-09-2017.fasta.gz -1 {file} --output outfile".format(
            sample_id=sample_id, file=file, outfile=outfile
        )
        out = subprocess.check_output(
            [
                "mykrobe",
                "genotype",
                sample_id,
                "data/tb-k21-probe-set-feb-09-2017.fasta.gz",
                "-1",
                file,
                "--output",
                outfile,
            ]
        )
        ## Load the output
        results = load_json(outfile)
        return results
