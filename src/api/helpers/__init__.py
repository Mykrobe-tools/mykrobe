import json


def load_json(f):
    with open(f, "r") as infile:
        return json.load(infile)
