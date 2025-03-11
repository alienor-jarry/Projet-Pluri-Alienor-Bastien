import math
import os
import sys
import numpy as np
from scipy.integrate import *
from matplotlib import pylab as plt
import bempp.api                    ## Need to add first, second, third and fourth derivate of BEMPP
import gmsh
import pyopencl 
import RK45
from einzel import einzel_lens


el0=einzel_lens(show_mesh=True)
el0.set_parameters(mesh_step=1.5, r1=1, d1=2,d2=2,t1=28,t2=26,t3=32)
el0.Ite(p=1,var=[27000,20,20])
