# singleton class Utility that has some utility functions
# to have always same form/colors in plots etc.:
# - create grid
# - TUGraz colors
# author: Markus Markl
# created: 07.11.2022
import subprocess


class Utility:
    """Singleton utility class for colors and adding grid lines to plots."""

    # colors of TU Graz presentation template
    col_tug = "#f70146"
    col_green = "#78b743"
    col_blue = "#285f82"
    col_yellow = "#e59352"
    col_cyan = "#77babf"
    col_purple = "#6c2f91"

    majorgridlw = 1.5
    minorgridlw = 1.0

    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Utility, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        pass

    def add_grid_to_axis(self, axis):
        axis.grid(which="major", color="#DDDDDD", linewidth=self.majorgridlw)
        axis.grid(which="minor", color="#EEEEEE", linewidth=self.minorgridlw, ls=":")
        axis.set_axisbelow(True)

    def get_git_version(self):
        try:
            git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode("utf-8")
            return git_hash
        except subprocess.CalledProcessError:
            return None


# Backwards compatibility alias
utility = Utility
