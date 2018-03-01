# Christopher Iliffe Sprague
# christopher.iliffe.sprague@gmail.com

import rov_dynamics
from dyn import System

class ROV(System):

    def __init__(self, tpr, b, r):
        System.__init__(self, [0, 0], [5, 5], 3, [tpr, b, r], rov_dynamics)
