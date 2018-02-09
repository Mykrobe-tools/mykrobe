import os


class Panel(object):

    def __init__(self, filepath, type="sequence"):
        self.filepath = filepath
        self.name = os.path.basename(filepath).split('.')[0]
        self.type = type

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name
