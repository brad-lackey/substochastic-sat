#!/usr/bin/python

from optimizeLUT import makeLUT, sendEmail,  parseTXT, LOOP_PENALTY
from subprocess32 import check_call, TimeoutExpired
import numpy as np
from joblib import Parallel, delayed
from itertools import product
from cleanupBrute import cleanup
import sys
import time

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: ./bruteOptimization.py <datfile>")
        sys.exit(1)

    # Use all CPUs minus 1
    N_JOBS = 2

    # Save at most this number of LUTs
    MAX_LUT = 10

    # Bins to divide the LUT
    bins = 5

    # Values of an A point in the LUT
    Avals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # list of A combinations
    A_list = [np.array(i) for i in product(Avals, repeat=bins)]

    dT = np.ones(bins)

    def saveProgress(progfile, index, success=True):
        with open(progfile, 'a') as f:
            if success:
                f.write("Job {0}/{1} Done\n".format(index + 1, len(A_list)))
            else:
                f.write("Job {0}/{1} Timed Out!\n".format(index+1, len(A_list)))

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
    indices = range(len(A_list))

    # delete finished indices from index list
    for i in sorted(done, reverse=True):
        del indices[i]

    class FileLock:
        def __init__(self):
            self.locked = False

        def lock(self):
            while self.locked:
                pass
            self.locked = True

        def unlock(self):
            self.locked = False

    # create the queue for results saving
    reslock = FileLock()

    # file to store the top MAX_LUT results
    resFile = tag + ".RESULTS.txt"

    CUTOFF_TIME = 60

    """Returns a dictionary of the MAX_LUT best loops"""
    def getResults(lock):
        results = {}

        lock.lock()

        try:
            with open(resFile, 'r') as f:
                line = f.readline()
                while(len(line) > 0):
                    iStr, loopStr, tStr, Astr = line.split(";")
                    i = int(iStr.split('=')[1])
                    loop = float(loopStr.split('=')[1])
                    t = float(tStr.rstrip('s').split('=')[1])

                    # save the loop from the file
                    results[i] = (loop, t)

                    line = f.readline()
        except IOError:
            pass  # No results yet

        lock.unlock()

        return results


    """Returns a dictionary of the MAX_LUT best loops"""
    def updateResults(index, loops, timeout, lock):
        results = {}

        lock.lock()

        new_min_found = False
        try:
            with open(resFile, 'r') as f:
                line = f.readline()
                while(len(line) > 0):
                    tup = line.split(";")
                    if len(tup) != 4:
                        line = f.readline()
                        continue

                    iStr, loopStr, tStr, Astr = tup
                    i = int(iStr.split('=')[1])
                    loop = float(loopStr.split('=')[1])
                    t = float(tStr.rstrip('s').split('=')[1])

                    if i == index:
                        line = f.readline()  # skip if i is given index
                        continue

                    # save the loop from the file
                    results[i] = (loop, t)

                    if not new_min_found and loops < loop:
                        new_min_found = True

                    line = f.readline()
        except IOError:
            pass  # No results yet

        if len(results) < MAX_LUT and index not in results:
            new_min_found = True

        results[index] = (loops, timeout)

        # only write to results if there is reason to
        if new_min_found:
            with open(resFile, 'w') as f:
                # iterate through results in sorted order by loop time, keeping only the top MAX_LUT results
                for i, tup in sorted(results.iteritems(), key=lambda x: x[1][0])[0:MAX_LUT]:
                    loop, t = tup
                    f.write("index={0}; loops={1}; time={2}s; A={3}\n".format(i, loop, t, A_list[i]))
            print("New minimum ({0}) added, at A={1}. See ".format(loops, A_list[index]) + resFile + " for details.")

        # set CUTOFF_TIME to the largest timeout in the top MAX_LUT results
        lastVal = sorted(results.values(), key=lambda x: x[1], reverse=True)[0][1]

        global CUTOFF_TIME
        if lastVal < CUTOFF_TIME:
            CUTOFF_TIME = lastVal
            print("Using new timeout of {0}".format(CUTOFF_TIME))

        lock.unlock()

        return results


    def bruteOptimize(index, A, queue):
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
        begin = time.time()
        try:
            check_call(args, timeout=CUTOFF_TIME)
        except TimeoutExpired:
            print("Job {0}/{1} Timed Out!".format(index+1, len(A_list)))
            saveProgress(progfile, index, success=False)  # save the index so that we don't have to redo it
            return LOOP_PENALTY
        timeout = time.time() - begin

        txtfile = fulltag + ".txt"
        hits, loops = parseTXT(txtfile)

        if hits < 1:
            loops = LOOP_PENALTY

        cleanup(fulltag, bins)

        saveProgress(progfile, index)

        updateResults(index, loops, timeout, queue)

        print("Job {0}/{1} Done".format(index+1, len(A_list)))

        return loops

    try:
        res = Parallel(n_jobs=N_JOBS, verbose=5)(delayed(bruteOptimize)(i, A_list[i], reslock) for i in indices)

        sendEmail("Optimization Finished!")

        # Print and save the best A's
        results = getResults(reslock)

        for rank, tup in enumerate(sorted(results.iteritems(), key=lambda x: x[1][0])):
            i = tup[0]
            loop, t = tup[1]
            print("Rank {0}: loops={1}, A={2}, time={3}".format(rank+1, loop, A_list[i], t))
            makeLUT(tag + ".LUT.{0}.BEST.{1}.txt".format(bins, rank+1), bins, dT, A_list[i])

    except KeyboardInterrupt:

        print("Cleaning output files...")

        # Cleans every output file up
        res = Parallel(n_jobs=-1)(delayed(cleanup)("{0}.{1}".format(tag, i), bins) for i, A in enumerate(A_list))
        print("Done!")

    sys.exit(0)
