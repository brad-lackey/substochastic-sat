#!/usr/bin/python

from optimizeLUT import tryLUT, makeLUT, sendEmail
import numpy as np
from joblib import Parallel, delayed
from itertools import product
import os
import sys

def cleanup(tag):
    try:
        os.remove(tag + ".LUT.txt")
    except Exception:
        pass  # No LUT file found
    try:
        os.remove(tag + ".log")
    except Exception:
        pass  # No log file found
    try:
        os.remove(tag + ".out")
    except Exception:
        pass  # No out file found
    try:
        os.remove(tag + ".txt")
    except Exception:
        pass  # No txt file found

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: ./cleanupBrute.py tag")
        sys.exit(1)

    tag = sys.argv[1]

    # Use all CPUs minus 1
    N_JOBS = -2

    # Bins to divide the LUT
    bins = 5

    # Values of an A point in the LUT
    Avals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # list of A combinations
    queue = [np.array(i) for i in product(Avals, repeat=bins)]

    dT = np.ones(bins)

    # Cleans every output file up
    res = Parallel(n_jobs=N_JOBS)(delayed(cleanup)("{0}.{1}".format(tag, i)) for i, A in enumerate(queue))

    sys.exit(0)
