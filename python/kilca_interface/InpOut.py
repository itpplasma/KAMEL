from read_in import *
from change_opts import *
from save_file import *

class InpOut:

    def __init__(self):
        pass

    def write(self, from_path, to_path):
        """Write Blueprint file."""

        raw = read_in(from_path)
        raw = change_opts(raw, self.ind, self.data(), self.sep)
        save_file(raw, to_path + '/' + self.BLUEPRINT)