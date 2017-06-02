#!/usr/bin/python

from optimizeLUT import makeLUT, sendEmail, parseTXT, UPDATE_PENALTY
from subprocess32 import check_call, TimeoutExpired
import numpy as np
from joblib import Parallel, delayed
from itertools import product
from cleanupBrute import cleanup
import sys
import time

if __name__ == "__main__":

    args = sys.argv
    verbosity = 0
    if '-v' in args:
        try:
            v = args[args.index('-v') + 1]
            verbosity = int(v)
            args.remove('-v')
            args.remove(v)
            if verbosity < 0 or verbosity > 2:
                raise ValueError
        except (IndexError, ValueError):
            print("Usage: ./bruteOptimization.py <datfile> [-v <verbosity (0: default no msgs, 1: update on minimum, 2: every job)>]")
            sys.exit(1)

    if len(args) != 2:
        print("Usage: ./bruteOptimization.py <datfile> [-v <verbosity (0: default no msgs, 1: update on minimum, 2: every job)>]")
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

    datfile = args[1]
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

    """Returns a dictionary of the MAX_LUT best updates"""
    def getResults(lock):
        results = {}

        lock.lock()

        try:
            with open(resFile, 'r') as f:
                line = f.readline()
                while(len(line) > 0):
                    iStr, updateStr, tStr, Astr = line.split(";")
                    i = int(iStr.split('=')[1])
                    update = float(updateStr.split('=')[1])
                    t = float(tStr.rstrip('s').split('=')[1])

                    # save the update from the file
                    results[i] = (update, t)

                    line = f.readline()
        except IOError:
            pass  # No results yet

        lock.unlock()

        return results


    """Returns a dictionary of the MAX_LUT best updates"""
    def updateResults(index, updates, timeout, lock):
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

                    iStr, updateStr, tStr, Astr = tup
                    i = int(iStr.split('=')[1])
                    update = float(updateStr.split('=')[1])
                    t = float(tStr.rstrip('s').split('=')[1])

                    if i == index:
                        line = f.readline()  # skip if i is given index
                        continue

                    # save the update from the file
                    results[i] = (update, t)

                    if not new_min_found and updates < update:
                        new_min_found = True

                    line = f.readline()
        except IOError:
            pass  # No results yet

        if len(results) < MAX_LUT and index not in results:
            new_min_found = True

        results[index] = (updates, timeout)

        # only write to results if there is reason to
        if new_min_found:
            with open(resFile, 'w') as f:
                # iterate through results in sorted order by loop time, keeping only the top MAX_LUT results
                for i, tup in sorted(results.iteritems(), key=lambda x: x[1][0])[0:MAX_LUT]:
                    update, t = tup
                    f.write("index={0}; updates={1}; time={2}s; A={3}\n".format(i, update, t, A_list[i]))

            if verbosity > 0:
                print("New minimum ({0}) added, at A={1}. See ".format(updates, A_list[index]) + resFile + " for details.")

        # set CUTOFF_TIME to the largest timeout in the top MAX_LUT results
        lastVal = sorted(results.values(), key=lambda x: x[1], reverse=True)[0][1]

        global CUTOFF_TIME
        if lastVal < CUTOFF_TIME:
            CUTOFF_TIME = lastVal

            if verbosity > 0:
                print("Using new timeout of {0}".format(CUTOFF_TIME))

        lock.unlock()

        return results


    def bruteOptimize(index, A, lock):
        fulltag = tag + "." + str(index)

        lut = fulltag + ".LUT.txt"

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
            if verbosity > 1:
                print("Job {0}/{1} Timed Out!".format(index+1, len(A_list)))
            saveProgress(progfile, index, success=False)  # save the index so that we don't have to redo it
            return UPDATE_PENALTY
        timeout = time.time() - begin

        txtfile = fulltag + ".txt"
        hits, updates = parseTXT(txtfile)

        if hits < 1:
            updates = UPDATE_PENALTY

        cleanup(fulltag, bins)

        saveProgress(progfile, index)

        updateResults(index, updates, timeout, lock)

        if verbosity > 1:
            print("Job {0}/{1} Done".format(index+1, len(A_list)))

        return updates

    try:
        res = Parallel(n_jobs=N_JOBS, verbose=5)(delayed(bruteOptimize)(i, A_list[i], reslock) for i in indices)

        # Print and save the best A's
        results = getResults(reslock)

        msg = "Optimization Finished!\n"

        for rank, tup in enumerate(sorted(results.iteritems(), key=lambda x: x[1][0])):
            i = tup[0]
            update, t = tup[1]
            msg += "Rank {0}: updates={1}, A={2}, time={3}\n".format(rank+1, update, A_list[i], t)
            makeLUT(tag + ".LUT.{0}.BEST.{1}.txt".format(bins, rank+1), bins, dT, A_list[i])

        print(msg)
        sendEmail(msg)

    except KeyboardInterrupt:

        print("Cleaning output files...")

        # Cleans every output file up
        res = Parallel(n_jobs=-1)(delayed(cleanup)("{0}.{1}".format(tag, i), bins) for i, A in enumerate(A_list))
        print("Done!")

    sys.exit(0)
