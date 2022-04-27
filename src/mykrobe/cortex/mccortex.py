from __future__ import print_function
import os
import sys
import subprocess
import logging
import tempfile

from mykrobe import K, MCCORTEX_BINARY_ENV_VAR

logger = logging.getLogger(__name__)


def syscall(command):
    if isinstance(command, list):
        command_str = " ".join(command)
    else:
        command_str = command

    logger.info(f"Run command: {command_str}")
    completed_process = subprocess.run(
        command,
        shell=not isinstance(command, list),
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logger.debug(f"Return code: {completed_process.returncode}")
    if completed_process.returncode != 0:
        print("Error running this command:", command_str, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise RuntimeError("Error in system call. Cannot continue")

    logger.debug(f"stdout:\n{completed_process.stdout.rstrip()}")
    logger.debug(f"stderr:\n{completed_process.stderr.rstrip()}")
    return completed_process


class McCortexRunner(object):
    def __init__(self):
        if MCCORTEX_BINARY_ENV_VAR in os.environ:
            self.mccortex31_path = os.environ[MCCORTEX_BINARY_ENV_VAR]
        else:
            dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
            self.mccortex31_path = os.path.join(dir_of_this_file, "mccortex31")

        if os.path.exists(self.mccortex31_path):
            pass
        elif os.path.exists(self.mccortex31_path + ".exe"):
            self.mccortex31_path += ".exe"
        else:
            raise RuntimeError(
                f"Did not find mccortex31. Expected it to be here: {self.mccortex31_path}. Cannot continue"
            )
        logger.debug(f"Found mccortex31: {self.mccortex31_path}")


class McCortexJoin(McCortexRunner):
    def __init__(
        self,
        sample,
        intersect_graph,
        ingraph,
    ):
        super(McCortexJoin, self).__init__()
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
        syscall(cmd)


class McCortexUnitigs(McCortexRunner):
    def __init__(self, ingraph):
        super(McCortexUnitigs, self).__init__()
        self.ingraph = ingraph

    def run(self):
        return self._run_cortex()

    def _run_cortex(self):
        cmd = [self.mccortex31_path, "unitigs", "-q", self.ingraph]
        return syscall(cmd)


class McCortexSubgraph(McCortexRunner):
    def __init__(
        self,
        sample,
        rmgraph,
        ingraph,
        tmp_dir=None,
    ):
        super(McCortexSubgraph, self).__init__()
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
        syscall(cmd)


class McCortexGenoRunner(McCortexRunner):
    def __init__(
        self,
        sample,
        panels,
        seq=None,
        ctx=None,
        kmer=K,
        threads=2,
        memory="1GB",
        force=False,
        panel_name=None,
        tmp_dir="tmp/",
        skeleton_dir="data/skeletons/",
    ):
        super(McCortexGenoRunner, self).__init__()
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
                "Force flag is false and coverage temp file exists"
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
            syscall(cmd)

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

            syscall(self.coverages_cmd)
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
