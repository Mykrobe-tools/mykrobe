from pkg_resources import get_distribution

try:
    __version__ = "v" + get_distribution("mykrobe").version
except:
    __version__ = "local"
