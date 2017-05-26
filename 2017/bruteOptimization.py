#!/usr/bin/python

from optimizeLUT import makeLUT, sendEmail,  parseTXT, LOOP_PENALTY
from subprocess32 import check_call, TimeoutExpired
import numpy as np
from joblib import Parallel, delayed
from itertools import product
from cleanupBrute import cleanup
import sys, os

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: ./bruteOptimization.py <datfile>")
        sys.exit(1)

    # Use all CPUs minus 1
    N_JOBS = 2

    # Bins to divide the LUT
    bins = 5

    # Values of an A point in the LUT
    Avals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # list of A combinations
    queue = [np.array(i) for i in product(Avals, repeat=bins)]

    dT = np.ones(bins)

    def saveProgress(progfile, index):
        with open(progfile, 'a') as f:
            f.write("Job {0}/{1} Done\n".format(index+1, len(queue)))

    def loadProgress(progfile):
        indices = []
        try:
            with open(progfile, 'r') as f:
                line = f.readline()
                while(len(line) > 0):
                    index = line.split()[1].split('/')[0]
                    indices.append(int(index)-1)  # save the 0-based index
                    line = f.readline()
        except IOError:
            pass  # No Progress file
        return indices


    datfile = sys.argv[1]
    # construct tag from datfile title
    tag = datfile.split('/')[-1].rstrip(".dat")
    progfile = tag + ".PROGRESS.txt"

    # get finished indices
    done = loadProgress(progfile)

    # create list of indices to complete
    indices = range(len(queue))

    # delete finished indices from index list
    for i in sorted(done, reverse=True):
        del indices[i]

    def bruteOptimize(index, A):
        fulltag = tag + "." + str(index)

        lut = fulltag + ".LUT." + str(bins) + ".txt"

        makeLUT(lut, bins, dT, A)

        args = []
        args.append('./testrun.pl')  # the program to run
        args.append('./ssmc')
        args.append(lut)
        args.append(datfile)
        args.append(str(1))
        args.append(fulltag)

        # returns 0 if successful otherwise throws error
        try:
            check_call(args, timeout=210)
        except TimeoutExpired:
            return LOOP_PENALTY

        txtfile = fulltag + ".txt"
        hits, loops = parseTXT(txtfile)

        if hits < 1:
            loops = LOOP_PENALTY

        cleanup(fulltag, bins)

        print("Job " + str(index+1) + "/" + str(len(queue)) + " Done")

        saveProgress(progfile, index)

        return loops

    try:
        res = Parallel(n_jobs=N_JOBS, verbose=5)(delayed(bruteOptimize)(i, queue[i]) for i in indices)

        sendEmail("Optimization Finished!")

        sortedRes = sorted(enumerate(res), key=lambda x:x[1])

        # Print and save the 10 best A's
        for i in range(10):
            index, loop = sortedRes[i]
            bestA = queue[index]
            print("Found " + str(loop) + ", A=" + str(bestA))
            makeLUT("a-h.2sat.LUT.BEST.{0}".format(i), 5, dT, bestA)

    except KeyboardInterrupt:

        print("Cleaning output files...")
        datfile = sys.argv[1]

        # construct tag from datfile title
        tag = datfile.split('/')[-1].rstrip(".dat")

        # Cleans every output file up
        res = Parallel(n_jobs=-1)(delayed(cleanup)("{0}.{1}".format(tag, i), bins) for i, A in enumerate(queue))
        print("Done!")

    sys.exit(0)
