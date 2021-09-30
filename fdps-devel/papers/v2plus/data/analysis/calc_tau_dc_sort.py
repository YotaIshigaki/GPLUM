import sys
from scipy.optimize import curve_fit
import numpy as np

def linear_fit(x, a):
    return x*a


n_smp = 100
np_1d = np.array([4, 8, 16, 32, 64])
np_2d = np_1d*np_1d
np_3d = np_2d*np_1d

time = np.array([3.53654e-05, 0.000124296, 0.0005277, 0.00188128, 0.00800101])

param, cov = curve_fit(linear_fit, np_2d*n_smp, time)
print(param, cov)
