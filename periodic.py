#Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize
#Lattice size


def periodic(z,coza,cozb):
    if z == 0:
        coza = 1
        cozb = Lz-1
    elif z == Lz-1:
        coza = 0
        cozb = Lz - 2
    else:
        coza = z + 1
        cozb = z - 1
    return(z,coza,cozb)
