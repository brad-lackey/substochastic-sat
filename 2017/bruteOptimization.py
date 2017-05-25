#!/usr/bin/python

from optimizeLUT import makeLUT, sendEmail,  parseTXT, LOOP_PENALTY
from subprocess32 import check_call, TimeoutExpired
import numpy as np
from joblib import Parallel, delayed
from itertools import product
from cleanupBrute import cleanup
import sys

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: ./bruteOptimization <datfile>")
        sys.exit(1)

    # Use all CPUs minus 1
    N_JOBS = 1

    # Bins to divide the LUT
    bins = 5

    # Values of an A point in the LUT
    Avals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # list of A combinations
    queue = [np.array(i) for i in product(Avals, repeat=bins)]

    dT = np.ones(bins)

    def bruteOptimize(index, A):
        tag = "a-h.2sat." + str(index)
        # loops = tryLUT(tag, "./ms_random/a-h.2sat.all.dat", 1, dT, A)

        lut = tag + ".LUT." + str(bins) + ".txt"

        makeLUT(lut, bins, dT, A)

        datfile = sys.argv[1]

        args = []
        args.append('./testrun.pl')  # the program to run
        args.append('./ssmc')
        args.append(lut)
        args.append(datfile)
        args.append(str(1))
        args.append(tag)

        # returns 0 if successful otherwise throws error
        try:
            check_call(args, timeout=1)
        except TimeoutExpired:
            return LOOP_PENALTY

        txtfile = tag + ".txt"
        hits, loops = parseTXT(txtfile)

        if hits < 1:
            loops = LOOP_PENALTY

        cleanup(tag, bins)

        print("Job " + str(index+1) + "/" + str(len(queue)) + " Done")
        # loops = index
        return loops

    res = Parallel(n_jobs=N_JOBS, verbose=5, backend="threading")(delayed(bruteOptimize)(i, A) for i, A in enumerate(queue))

    sendEmail("Optimization Finished!")

    sortedRes = sorted(enumerate(res), key=lambda x:x[1])

    # Print and save the 10 best A's
    for i in range(10):
        index, loop = sortedRes[i]
        bestA = queue[index]
        print("Found " + str(loop) + ", A=" + str(bestA))
        makeLUT("a-h.2sat.LUT.BEST.{0}".format(i), 5, dT, bestA)

    # Cleans every output file up
    res = Parallel(n_jobs=N_JOBS)(delayed(cleanup)("a-h.2sat.{0}".format(i), bins) for i, A in enumerate(queue))

    sys.exit(0)