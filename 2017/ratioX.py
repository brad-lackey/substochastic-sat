#!/usr/bin/python

from optimizeLUT import parseLUT, sendEmail, plotLUT, plotPsize, tryLUT
from histAnalysis import parseOUT
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print "Usage: ./ratioX.py <LUT> <DAT> tag"
        sys.exit(1)

    lut = sys.argv[1]
    datfile = sys.argv[2]
    tag = sys.argv[3]

    optima = {}
    ratio = {}

    with open(datfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            fname, oLbl, _, oStr, tLbl, _, tStr = line.split()

            optimum = oStr.split('=')[-1].strip()
            optima[fname] = optimum

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

    tryLUT(tag, datfile, 1, dT, A, psize)

    outfile = tag + ".out"

    files, times, loops, updates = parseOUT(outfile)

    x = []
    y = []
    for i in range(len(files)):
        x.append(ratio[files[i]])
        y.append(updates[i])

    print("Ratio:{0}\nUpdates:{1}\n".format(x, y))

    plt.plot(x, y, 'bo')
    plt.xlabel("Ratio")
    plt.ylabel("Updates")
    plt.title("Updates vs. Var/Clause Ratio")
    plt.show()