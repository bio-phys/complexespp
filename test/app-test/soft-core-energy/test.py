import sys
import pandas as pd
import numpy as np

stat = pd.read_csv('test.stat')
ref = np.loadtxt('reference.dat')

res = True
try:
    np.testing.assert_almost_equal(ref, stat['softcore-pot']+stat['softcore-el'], decimal=4)
except:
    res = False

print("all ok? = {}".format(res))
sys.exit(not res)