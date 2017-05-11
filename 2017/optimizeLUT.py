#!/usr/bin/python
import sys
from createLUT import makeLUT
from subprocess import check_call
import numpy as np

"""Returns the percentage of hits and avg runtime as a tuple"""
def parseTXT(txtfile):
    last = ''
    with open(txtfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            last = line
            line = f.readline()

    # Parse last line
    _, hitStr, _, tStr = last.split()

    fraction = map(float, hitStr[0:hitStr.rfind('(')].split('/'))
    hit = fraction[0]/fraction[1]

    t = float(tStr.rstrip("s"))

    return hit, t

"""Returns the avg time of a set of conf files using given LUT"""
def tryLUT(tag, filename, trials, dT, A, weight, runtime):
    if len(dT) != len(A):
        raise Exception("Vectors dT and A are not the same length!")

    lut = tag + ".LUT.txt"
    bins = len(dT)

    makeLUT(lut, bins, dT, A)

    args = []
    args.append('./testrun.pl') # the program to run
    args.append('./ssmc')
    args.append(lut)
    args.append(filename)
    args.append(trials)
    args.append(tag)

    if weight and runtime:
        args.append(weight)
        args.append(runtime)

    # returns 0 if successful otherwise throws error
    check_call(args)

    txtfile = tag + ".txt"
    hits, avg_time = parseTXT(txtfile)

    penalty = 10000
    if hits < 1:
        avg_time = penalty

    return avg_time


def main():
    if len(sys.argv) == 5 or len(sys.argv) == 7:
        filename = sys.argv[1]
        bins = int(sys.argv[2])
        trials = sys.argv[3]
        tag = sys.argv[4]
        weight = None
        runtime = None

        if len(sys.argv) == 7:
            weight = sys.argv[5]
            runtime = sys.argv[6]

    else:
        print("Usage: ./optimizeLUT <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
        return 1

    # Guesses for dT and A
    dT = np.ones(bins)
    A = np.linspace(1, 0, bins)

    # Just a random test
    print(tryLUT(tag, filename, trials, dT, A, weight, runtime))

    return 0

if __name__ == "__main__":
    sys.exit(main())