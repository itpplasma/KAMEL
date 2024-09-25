from .read_in import *
from .change_opts import *
from .save_file import *

class InpOut:

    def __init__(self):
        pass

    def write(self, from_path, to_path):
        """Write Blueprint file."""

        raw = read_in(from_path)
        raw = change_opts(raw, self.ind, self.return_data(), self.sep)
        save_file(raw, to_path + '/' + self.BLUEPRINT)

    
    def return_data(self):
        """
        Description:
            Return the data of the class used to write the corresponding .in file.
        """
        l = []
        for key in self.data:
            l.append(self.data[key])
        return l
