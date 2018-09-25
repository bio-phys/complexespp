import sys
import pandas as pd
import numpy as np

stat = pd.read_csv('test.stat')
ref = np.loadtxt('reference.dat')

res = True
try:
    np.testing.assert_almost_equal(ref, stat.WCA, decimal=5)
except:
    res = False

print("all ok? = {}".format(res))
sys.exit(not res)