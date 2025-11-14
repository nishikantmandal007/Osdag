import numpy as np
import math

# This file is intentionally simple.
# The complex, custom PSO algorithm was removed.
# We now use the standard 'pyswarm' library's 'pso' function,
# which is called directly from 'weldedPlateGirder.py'.

class Section:
    """Helper class to store section properties for optimization particles."""
    def __init__(self):
        self.tf = self.tw = self.bf = self.D = None
        self.tf_top = self.tf_bot = self.bf_top = self.bf_bot = None
        self.c = self.t_stiff = None