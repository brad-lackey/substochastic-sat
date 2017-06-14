#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
from subprocess32 import check_call
import time
from optimizeLUT import sendEmail

def parseOUT(filename):

    files = []
    times = []
    loops = []
    updates = []

    with open(filename, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            file, _, tStr, loopStr, uStr = line.split()

            files.append(file)
            times.append(float(tStr))
            loops.append(int(loopStr))
            updates.append(int(uStr))

            line = f.readline()

    return files, times, loops, updates

if __name__ == "__main__":

    if len(sys.argv) != 5:
        print "Usage: ./histAnalysis.py <LUT> <DAT> trials tag"
        sys.exit(1)

    lut = sys.argv[1]
    datfile = sys.argv[2]
    trials = int(sys.argv[3])
    tag = sys.argv[4]

    args = []
    args.append('./testrun.pl')  # the program to run
    args.append('./ssmc')
    args.append(lut)
    args.append(datfile)
    args.append(str(trials))
    args.append(tag)

    begin = time.time()
    # run testrun.pl with the given args
    check_call(args)
    elapsed = time.time() - begin

    print("Elapsed: {0} seconds".format(elapsed))

    sendEmail("Histogram analysis done.")

    _, times, loops, updates = parseOUT(tag + ".out")

    plt.figure(0)
    plt.hist(times, 100)
    plt.title("Times Histogram for {0}".format(tag))

    plt.figure(1)
    plt.hist(loops, 100)
    plt.title("Loops Histogram for {0}".format(tag))

    plt.figure(2)
    plt.hist(updates, 100)
    plt.title("Updates Histogram for {0}".format(tag))

    plt.show()

    sys.exit(0)