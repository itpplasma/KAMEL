# singleton class utility that has some utility functions 
# to have always same form/colors in plots etc.:
# - create grid
# - TUGraz colors
# author: Markus Markl
# created: 07.11.2022

class utility:
    # colors of TU Graz presentation template
    col_tug = '#f70146'
    col_green = '#78b743'
    col_blue = '#285f82'
    col_yellow = '#e59352'
    col_cyan = '#77babf'
    col_purple = '#6c2f91'

    majorgridlw = 1.5
    minorgridlw = 1.0

    def __init__(self):
        pass

    def __call__(self, cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super(utility, cls).__call__(cls)
        return cls.instance

    def add_grid_to_axis(self, axis):
        axis.grid(which='major', color='#DDDDDD', linewidth=self.majorgridlw)
        axis.grid(which='minor', color='#EEEEEE',
                  linewidth=self.minorgridlw, ls=':')
        axis.set_axisbelow(True)
