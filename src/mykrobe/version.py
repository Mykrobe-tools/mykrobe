from importlib.metadata import PackageNotFoundError
from importlib.metadata import version


try:
    __version__ = f"v{version('mykrobe')}"
except PackageNotFoundError:
    __version__ = "local"
