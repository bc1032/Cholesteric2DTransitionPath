#Libraries
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import scipy
from scipy.optimize import minimize

import periodic
import calculategrad
import initialise
import update
import energy as en
import analyse
import writeenergy

def secondvariation(guess,original,GradE,sig,Lx,Lz,Lt,ks,kt,q0,z,t,s,alpha,beta,gamma,a,b,c,Q1,Q2,Q3,Q4,Q5,ws,timeen,splay,twist,bend,surface,bulk):
