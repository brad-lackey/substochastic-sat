#!/usr/bin/python

from optimizeLUT import tryLUT
from utilities import parseLUT, parseOUT
import matplotlib.pyplot as plt
import sys
import numpy as np

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print "Usage: ./ratioX.py <LUT> <DAT> tag"
        sys.exit(1)

    lut = sys.argv[1]
    datfile = sys.argv[2]
    tag = sys.argv[3]

    optima = {}
    times = {}
    ratio = {}

    with open(datfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            fname, oLbl, _, oStr, tLbl, _, tStr = line.split()

            optimum = int(oStr.split('=')[-1].strip())
            optima[fname] = optimum
            t = float(tStr.split('=')[-1].strip())
            times[fname] = t


            line = f.readline()

    for cnf in optima.keys():
        with open(cnf, 'r') as f:
            line = f.readline()
            while len(line) > 0:
                if line.startswith('p'):
                    p, fmt, varStr, cStr = line.split()

                    var = float(varStr)
                    clauses = float(cStr)

                    ratio[cnf] = var/clauses

                line = f.readline()


    bins, dT, A, psize = parseLUT(lut)

    tryLUT(var, tag, datfile, 1, dT, A, psize)

    outfile = tag + ".out"

    files, _optima, _times, _loops, updates = parseOUT(outfile)

    x = []
    y = []
    for i in range(len(files)):
        x.append(ratio[files[i]])
        y.append(updates[i])

    x = np.array(x)

    hist, edges = np.histogram(x, len(files))

    avgs = (edges[1:] + edges[:-1])/2

    for i in range(len(hist)):
        if hist[i] > 0:
            if i == len(hist)-1:
                indices = np.flatnonzero(np.logical_and(x >= edges[i], x <= edges[i+1]))
            else:
                indices = np.flatnonzero(np.logical_and(x >= edges[i], x < edges[i+1]))
            datfile = tag + ".{0:.3f}.dat".format(avgs[i])

            with open(datfile, 'w') as f:
                for j in indices:
                    f.write(files[j])
                    f.write("   ")
                    f.write("O = {0}   T = {1}\n".format(optima[files[j]], times[files[j]]))





    print("Ratio:{0}\nUpdates:{1}\n".format(x, y))

    plt.plot(x, y, 'bo')
    plt.xlabel("Ratio")
    plt.ylabel("Updates")
    plt.title("Updates vs. Var/Clause Ratio")
    plt.show()
