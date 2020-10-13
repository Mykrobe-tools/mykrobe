import json
import logging
import os
import re
import requests
import shutil
import tarfile

from mykrobe.species_data import SpeciesDir

logger = logging.getLogger(__name__)

MANIFEST_URL = "https://raw.githubusercontent.com/Mykrobe-tools/mykrobe-data/main/mykrobe_panels_manifest.json"

class DataDir:
    def __init__(self, root_dir):
        self.root_dir = os.path.abspath(root_dir)
        self.manifest_json = os.path.join(self.root_dir, "manifest.json")
        self.lock_file = os.path.join(self.root_dir, ".lock")
        self.load_manifest()

    def is_locked(self):
        return os.path.exists(self.lock_file)

    def start_lock(self, error_message=None):
        if self.is_locked():
            if error_message is not None:
                raise RuntimeError(error_message)
            else:
                raise RuntimeError(f"Lock file found: {self.lock_file}")
        with open(self.lock_file, "w"):
            pass

    def stop_lock(self):
        assert self.is_locked()
        os.unlink(self.lock_file)

    def load_manifest(self):
        if os.path.exists(self.manifest_json):
            with open(self.manifest_json) as f:
                self.manifest = json.load(f)
        else:
            self.manifest = {}

    def create_root(self):
        if not os.path.exists(self.root_dir):
            logger.info(f"Creating panels directory {self.root_dir}")
            os.mkdir(self.root_dir)

    def save_manifest(self):
        self.create_root()
        with open(self.manifest_json, "w") as f:
            json.dump(self.manifest, f, indent=2, sort_keys=True)

    def update_manifest(self, filename=None, url=None):
        self.create_root()
        error_message = f"Cannot update panel list. Looks like another process is already trying to modify the mykrobe data directory {self.root_dir}. If this is not the case, then you can delete this file and try again: {self.lock_file}"
        self.start_lock(error_message)
        if filename is not None:
            logger.info(f"Getting latest panel information from file {filename}")
            with open(filename) as f:
                new_manifest = json.load(f)
        else:
            if url is None:
                url = MANIFEST_URL
            try:
                logger.info(f"Getting panels information from {url}")
                new_manifest = json.loads(requests.get(url).text)
            except:
                raise RuntimeError(f"Error getting latest panel information from {url}")

        for species, species_dict in new_manifest.items():
            logger.info(f"Updating metadata for species {species}")
            if species not in self.manifest:
                self.manifest[species] = {"installed": None}
            self.manifest[species]["latest"] = species_dict

        self.save_manifest()
        self.stop_lock()
        logger.info(f"Finished updating metadata in panels directory {self.root_dir}")

    def print_panels_summary(self):
        if len(self.manifest) == 0:
            print("No data")
            return

        species_lines = []
        print("\nSpecies summary:\n")
        print("Species", "Update_available", "Installed_version", "Installed_url", "Latest_version", "Latest_url", sep="\t")
        for species, data in sorted(self.manifest.items()):
            if self.species_is_up_to_date(species):
                update_available = "no"
            else:
                update_available = "yes"

            if data["installed"] is None:
                installed_version = "None"
                installed_url = "NA"
            else:
                installed_version = data["installed"]["version"]
                installed_url = data["installed"]["url"]

            if data["installed"] is not None:
                sdir = self.get_species_dir(species)
                species_lines.append(f"\n{species} default panel: {sdir.default_panel()}")
                species_lines.append(f"{species} panels:")
                species_lines.append(f"Panel\tReference\tDescription")
                for panel in sdir.panel_names():
                    sdir.set_panel(panel)
                    species_lines.append("\t".join([panel, sdir.reference_genome(), sdir.description()]))

            print(species, update_available, installed_version, installed_url, data["latest"]["version"], data["latest"]["url"], sep="\t")

        if len(species_lines):
            print(*species_lines, sep="\n")
        else:
            print("\nNo panels are installed")


    def add_or_replace_species_data(self, tarball_name, force=False):
        error_message = f"Cannot add/update new species. Looks like another process is already trying to modify the mykrobe data directory {self.root_dir}. If this is not the case, then you can delete this file and try again: {self.lock_file}"
        self.start_lock(error_message)
        tmp_dir = os.path.join(self.root_dir, "tmp.add_species")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)

        # For testing, is useful to install from files. But this is hidden from
        # the user because we wan't them to download out tarballs. So check
        # if the tarball_name looks like a url, and if not assume it's a file.
        # Note: this regex is not generally correct for deciding a url or
        # not, but our downloads will always match this, so is ok here.
        from_file = re.search(r"^https?://", tarball_name) is None
        if from_file:
            logger.info(f"Installing species from file {tarball_name}")
            to_extract = tarball_name
        else:
            logger.info(f"Downloading file {tarball_name}")
            request = requests.get(tarball_name, allow_redirects=True)
            to_extract = os.path.join(tmp_dir, "tmp.add_species.tar.gz")
            with open(to_extract, 'wb') as t:
                t.write(request.content)

        logger.info(f"Extracting tarball {to_extract}")
        with tarfile.open(to_extract, mode="r") as t:
            t.extractall(path=tmp_dir)

        if not from_file:
            os.unlink(to_extract)

        dirs = os.listdir(tmp_dir)
        assert len(dirs) == 1
        extracted_dir = os.path.join(tmp_dir, dirs[0])
        spdir = SpeciesDir(extracted_dir)
        if not spdir.sanity_check():
            raise RuntimeError(f"Something wrong with data in {extracted_dir}. Cannot continue")

        species = spdir.species_name()
        if self.species_is_installed(species) and not force:
            raise RuntimeError(f"Species '{species}', already exists. To replace with new version, use the force option. Cannot continue")

        new_dir = os.path.join(self.root_dir, species)
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
        os.rename(extracted_dir, new_dir)
        shutil.rmtree(tmp_dir)
        self.manifest[species]["installed"] = {"version": spdir.version(), "url": tarball_name}
        self.save_manifest()
        self.stop_lock()



    def remove_species(self, species):
        if not self.species_is_installed(species):
            raise ValueError(f"Species {species} is not installed, so cannot remove it")
        error_message = f"Cannot remove species. Looks like another process is already trying to modify the mykrobe data directory {self.root_dir}. If this is not the case, then you can delete this file and try again: {self.lock_file}"
        self.start_lock(error_message)
        self.manifest[species]["installed"] = None
        species_directory = os.path.join(self.root_dir, species)
        if os.path.exists(species_directory):
            shutil.rmtree(species_directory)
        self.save_manifest()
        self.stop_lock()

    def all_species_list(self):
        return sorted(list(self.manifest.keys()))

    def installed_species(self):
        return [x for x in self.all_species_list() if self.species_is_installed(x)]

    def species_is_installed(self, species):
        return species in self.manifest and self.manifest[species]["installed"] is not None

    def species_is_up_to_date(self, species):
        if not self.species_is_installed(species):
            return False

        d = self.manifest[species]
        return d["installed"] is not None and int(d["latest"]["version"]) <= int(d["installed"]["version"])

    def get_species_dir(self, species_name):
        if species_name not in self.manifest:
            raise ValueError(f"Species '{species_name}' not found in data directory {self.root_dir}")
        elif not self.species_is_installed(species_name):
            return None
        else:
            return SpeciesDir(os.path.join(self.root_dir, species_name))

    def update_species(self, species):
        if species not in self.manifest:
            raise KeyError(f"Unknown species '{species}'. Must choose from: {self.all_species_list()}")
        if self.species_is_up_to_date(species):
            logger.info(f"Species {species} is up to date. Nothing to do.")
        else:
            tarball = self.manifest[species]["latest"]["url"]
            logger.info(f"{species}: updating from {tarball}")
            self.add_or_replace_species_data(tarball, force=True)
            logger.info(f"Updated species {species}")

    def update_all_species(self):
        logger.info("Updating all species")
        for species in self.all_species_list():
            self.update_species(species)
        logger.info("Finished updating")
