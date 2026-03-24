import importlib.metadata as metadata
import importlib.util
from pathlib import Path


VERSION_MODULE = Path(__file__).resolve().parents[1] / "src" / "mykrobe" / "version.py"


def load_version_module():
    spec = importlib.util.spec_from_file_location("mykrobe_version_test", VERSION_MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_version_uses_installed_package_metadata(monkeypatch):
    monkeypatch.setattr(metadata, "version", lambda package_name: "1.2.3")

    module = load_version_module()

    assert module.__version__ == "v1.2.3"


def test_version_falls_back_to_local_when_package_not_installed(monkeypatch):
    def raise_package_not_found(package_name):
        raise metadata.PackageNotFoundError

    monkeypatch.setattr(metadata, "version", raise_package_not_found)

    module = load_version_module()

    assert module.__version__ == "local"
