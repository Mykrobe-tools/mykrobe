import json
import logging
import os

logger = logging.getLogger(__name__)

class SpeciesDir:
    def __init__(self, root_dir):
        self.root_dir = os.path.abspath(root_dir)
        if not os.path.exists(self.root_dir):
            raise FileNotFoundError(f"Species directory {self.root_dir} not found.")
        self.manifest_json = os.path.join(root_dir, "manifest.json")
        if not os.path.exists(self.manifest_json):
            raise FileNotFoundError(f"Manifest file not found in species directory {self.root_dir}. Expected to find {self.manifest_json}.")
        self.panel = None
        with open(self.manifest_json) as f:
            self.manifest = json.load(f)
        self.set_panel(self.default_panel())

    def set_panel(self, panel_name):
        if panel_name not in self.manifest["panels"]:
            raise KeyError(f"Panel '{panel_name}' not found from species directory {self.root_dir}. Available panels: {self.panel_names()}")
        self.panel_name = panel_name
        self.panel = self.manifest["panels"][self.panel_name]

    def default_panel(self):
        return self.manifest["default_panel"]

    @classmethod
    def _check_keys_in_dict(cls, keys_to_check, d):
        all_ok = True
        for key in keys_to_check:
            key_found = key in d
            logger.debug(f"Check '{key}' found: {key_found}")
            if not key_found:
                all_ok = False
        return all_ok

    def sanity_check(self):
        logger.debug(f"Checking all panels in {self.root_dir}")
        expect_keys = ["species_name", "version", "panels", "default_panel"]
        all_ok = SpeciesDir._check_keys_in_dict(expect_keys, self.manifest)
        if "panels" not in self.manifest:
            logger.warning(f"No 'panels' section found in JSON file {self.manifest_json}")
            return False
        expect_keys = ["description", "reference_genome", "fasta_files", "json_files", "kmer", "species_phylo_group"]
        json_types = ["amr", "lineage", "hierarchy"]

        for panel_name in self.panel_names():
            logger.debug(f"Checking panel {panel_name} in {self.root_dir}")
            self.set_panel(panel_name)
            all_ok = all_ok and SpeciesDir._check_keys_in_dict(expect_keys, self.panel)
            fasta_files = self.fasta_files()
            if fasta_files is None or len(fasta_files) == 0:
                logger.warning(f"Must have at least one FASTA file name in {self.manifest_json}")
                all_ok = False
            else:
                for filename in fasta_files:
                    file_found = os.path.exists(filename)
                    logger.debug(f"Found file {filename}: {file_found}")
                    all_ok = all_ok and file_found
                    if not file_found:
                        logger.warning(f"FASTA file not found: {filename}")

            found_at_least_one_json = False
            for json_type in json_types:
                json_file = self.json_file(json_type)
                if json_file is None:
                    continue
                if  os.path.exists(json_file):
                    found_at_least_one_json = True
                else:
                    logger.warning(f"JSON {json_type} file not found: {json_file}")

            all_ok = all_ok and found_at_least_one_json

        return all_ok

    def panel_names(self):
        return sorted(list(self.manifest["panels"].keys()))

    def species_name(self):
        return self.manifest["species_name"]

    def version(self):
        return self.manifest["version"]

    def description(self):
        return self.panel["description"]

    def reference_genome(self):
        return self.panel["reference_genome"]

    def kmer(self):
        return self.panel["kmer"]

    def species_phylo_group(self):
        return self.panel["species_phylo_group"]

    def _absolute_path(self, filename):
        if filename is None:
            return None
        else:
            return os.path.join(self.root_dir, filename)

    def fasta_files(self):
        if self.panel.get("fasta_files", None) is None:
            return None
        else:
            return [self._absolute_path(f) for f in self.panel["fasta_files"]]

    def json_file(self, data_type):
        try:
            filename = self.panel["json_files"][data_type]
        except:
            return None
        if filename is None:
            return None
        else:
            return self._absolute_path(filename)

