from __future__ import print_function
import os
import sys
import subprocess
import logging
import tempfile

logger = logging.getLogger(__name__)


class McCortexRunner(object):
    def __init__(self, mccortex31_path):
        self.mccortex31_path = mccortex31_path

    @property
    def mccortex31_path(self):
        if not hasattr(self, "_mccortex31_path"):
            logger.error("mccortex31_path attribute used before setter called")
            exit(-1)
        return self._mccortex31_path

    @mccortex31_path.setter
    def mccortex31_path(self, suggested_fpath):
        logger.debug("Setting path to mccortex executable")
        if os.path.isfile(suggested_fpath):
            logger.debug(
                "Suggested path to mccortex executable is valid: %s", suggested_fpath
            )
            self._mccortex31_path = suggested_fpath
            return

        logger.debug(
            "Suggested path to mccortex executable is invalid: %s", suggested_fpath
        )
        fallback_path = os.path.join(os.path.dirname(sys.argv[0]), "mccortex31")
        logger.debug(
            "Checking fallback mccortex executable file path: %s", fallback_path
        )

        if not os.path.isfile(fallback_path):
            logger.debug(
                "Fallback mccortex executable file path invalid: %s", fallback_path
            )
            exit(-1)

        logger.debug("Valid path to mccortex executable found: %s", fallback_path)
        self._mccortex31_path = fallback_path


class McCortexJoin(McCortexRunner):
    def __init__(self, sample, intersect_graph, ingraph, mccortex31_path="mccortex31"):
        super(McCortexJoin, self).__init__(mccortex31_path)
        self.sample = sample
        self.intersect_graph = intersect_graph
        self.ingraph = ingraph
        self.out_ctx_dir = tempfile.mkdtemp()
        self.out_ctx_path = os.path.join(self.out_ctx_dir, "%s_int.ctx" % sample)

    def run(self):
        self._run_cortex()
        return self.out_ctx_path

    def _run_cortex(self):
        cmd = [
            self.mccortex31_path,
            "join",
            "-q",
            "--out",
            self.out_ctx_path,
            "--intersect",
            self.intersect_graph,
            self.ingraph,
        ]
        subprocess.check_output(cmd)


class McCortexUnitigs(McCortexRunner):
    def __init__(self, ingraph, mccortex31_path="mccortex31"):
        super(McCortexUnitigs, self).__init__(mccortex31_path)
        self.ingraph = ingraph

    def run(self):
        return self._run_cortex()

    def _run_cortex(self):
        cmd = [self.mccortex31_path, "unitigs", "-q", self.ingraph]
        return subprocess.check_output(cmd)


class McCortexSubgraph(McCortexRunner):
    def __init__(
        self, sample, rmgraph, ingraph, tmp_dir=None, mccortex31_path="mccortex31"
    ):
        super(McCortexSubgraph, self).__init__(mccortex31_path=mccortex31_path)
        self.rmgraph = rmgraph
        self.sample = sample
        self.ingraph = ingraph
        if tmp_dir is None:
            self.out_ctx_dir = tempfile.mkdtemp()
        else:
            self.out_ctx_dir = tmp_dir
        self.out_ctx_path = os.path.join(self.out_ctx_dir, "%s_new.ctx" % sample)

    def run(self):
        self._run_cortex()
        return self.out_ctx_path

    def _run_cortex(self):
        cmd = [
            self.mccortex31_path,
            "view",
            "-q",
            "-k",
            self.rmgraph,
            "|",
            "awk",
            "'{print $1}'",
            "|",
            self.mccortex31_path,
            "subgraph",
            "-q",
            "--out",
            self.out_ctx_path,
            "--invert",
            "--seq",
            "-",
            self.ingraph,
        ]
        logger.debug(subprocess.list2cmdline(cmd))
        logger.debug(" ".join(cmd))
        cmd = subprocess.Popen(" ".join(cmd), shell=True).wait()


class McCortexGenoRunner(McCortexRunner):
    def __init__(
        self,
        sample,
        panels,
        seq=None,
        ctx=None,
        kmer=31,
        threads=2,
        memory="1GB",
        force=False,
        panel_name=None,
        tmp_dir="tmp/",
        skeleton_dir="data/skeletons/",
        mccortex31_path="mccortex31",
    ):
        super(McCortexGenoRunner, self).__init__(mccortex31_path)
        self.sample = sample
        self.panels = panels
        self.seq = seq
        self.ctx = ctx
        self.kmer = kmer
        self.force = force
        self._panel_name = panel_name
        self.tmp_dir = tmp_dir
        self.threads = threads
        self.memory = memory
        self.skeleton_dir = skeleton_dir
        if self.seq and self.ctx:
            raise ValueError("Can't have both -1 and -c")

    def run(self):
        if self.force or not os.path.exists(self.covg_tmp_file_path):
            self._check_panels()
            self._run_cortex()
        else:
            logger.warning(
                "Not running mccortex. "
                "Force flag is false or coverage temp file exists"
            )

    def _check_panels(self):
        # If panel does not exists then build it
        for panel in self.panels:
            if not os.path.exists(panel.filepath):
                raise ValueError("Could not find a panel at %s." % panel.filepath)

    def _run_cortex(self):
        # If ctx binary does not exist then build it
        self._build_panel_binary_if_required()
        # Now get coverage on panel
        self._run_coverage_if_required()

    def _build_panel_binary_if_required(self):
        if not os.path.exists(self.ctx_skeleton_filepath) or self.force:
            if os.path.exists(self.ctx_skeleton_filepath):
                os.remove(self.ctx_skeleton_filepath)
            # panel
            seq_list = self._create_sequence_list()
            cmd = (
                [
                    self.mccortex31_path,
                    "build",
                    "-q",
                    "-m %s" % self.memory,
                    "-t",
                    "%i" % self.threads,
                    "-k",
                    str(self.kmer),
                ]
                + seq_list
                + [self.ctx_skeleton_filepath]
            )
            logger.debug("Executing command:\n%s", cmd)
            subprocess.check_output(cmd)

    def _create_sequence_list(self):
        seq_list = []
        seq_list.extend(["-s", "%s" % self.panel_name[:100]])
        for panel in self.panels:
            seq_list.extend(["-1", panel.filepath])
        return seq_list

    @staticmethod
    def _execute_command(command):
        process = subprocess.Popen(command, stdout=subprocess.PIPE)

        while True:
            nextline = process.stdout.readline()
            if not nextline and process.poll() is not None:
                break
            sys.stdout.write(nextline.decode("utf-8"))

        output = process.communicate()[0]
        exit_code = process.returncode
        if exit_code == 0:
            return output
        else:
            raise subprocess.CalledProcessError

    def _run_coverage_if_required(self):
        if (
            not os.path.exists(self.ctx_tmp_filepath)
            or not os.path.exists(self.covg_tmp_file_path)
            or self.force
        ):
            if os.path.exists(self.ctx_tmp_filepath):
                os.remove(self.ctx_tmp_filepath)
            if os.path.exists(self.covg_tmp_file_path):
                os.remove(self.covg_tmp_file_path)

            logger.debug("Running coverages command:\n%s", " ".join(self.coverages_cmd))
            try:
                self._execute_command(self.coverages_cmd)
            except subprocess.CalledProcessError:
                command = subprocess.list2cmdline(self.coverages_cmd)
                exception_message = (
                    """mccortex31 raised an error. 
                    Is it on PATH? check by running `mccortex31 geno`. 
                    'The command that through the error was `%s` """
                    % command
                )
                raise ValueError(exception_message)
        else:
            logger.warning(
                "Using pre-built binaries. Run with --force if "
                "panel has been updated."
            )

    @property
    def coverages_cmd(self):
        if self.seq:
            return self.coverages_cmd_seq
        elif self.ctx:
            return self.coverages_cmd_ctx
        else:
            raise ValueError("Need either seq or ctx binary to run coverages")

    @property
    def base_geno_command(self):
        return [
            self.mccortex31_path,
            "geno",
            "-t",
            "%i" % self.threads,
            "-m %s" % self.memory,
            "-k",
            str(self.kmer),
            "-o",
            self.covg_tmp_file_path,
        ]

    @property
    def coverages_cmd_seq(self):
        cmd = self.base_geno_command
        cmd.extend(["-I", self.ctx_skeleton_filepath])
        cmd.extend(["-s", self.sample_name])
        for seq in self.seq:
            cmd.extend(["-1", seq])
        for panel in self.panels:
            cmd.extend(["-c", panel.filepath])
        cmd.append(self.ctx_tmp_filepath)
        return cmd

    @property
    def coverages_cmd_ctx(self):
        cmd = self.base_geno_command
        cmd.extend(["-g", self.ctx])
        for panel in self.panels:
            cmd.extend(["-c", panel.filepath])
        cmd.append(self.ctx_tmp_filepath)
        return cmd

    @property
    def sample_name(self):
        return "-".join([self.sample, str(self.kmer)])

    @property
    def panel_name(self):
        if self._panel_name is None:
            self._panel_name = "-".join([p.name.replace("/", "-") for p in self.panels])
        return self._panel_name

    @property
    def sample_panel_name(self):
        return "_".join([self.sample_name, self.panel_name])

    @property
    def ctx_tmp_filepath(self):
        sample_panel_name = self.sample_panel_name + ".ctx"
        return os.path.join(self.tmp_dir, sample_panel_name)

    @property
    def covg_tmp_file_path(self):
        sample_panel_name = self.sample_panel_name + ".covgs"
        return os.path.join(self.tmp_dir, sample_panel_name)

    @property
    def ctx_skeleton_filepath(self):
        panel_name = self.panel_name.replace("/", "-")[:100]
        panel_fname = "{panel_name}_{kmer}.ctx".format(
            panel_name=panel_name, kmer=self.kmer
        )
        return os.path.join(self.skeleton_dir, panel_fname)

    def remove_temporary_files(self):
        os.remove(self.ctx_tmp_filepath)
        os.remove(self.covg_tmp_file_path)
