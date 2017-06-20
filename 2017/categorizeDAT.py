#!/usr/bin/python

from optimizeLUT import parseLUT, sendEmail, plotLUT, plotPsize, tryLUT
import matplotlib.pyplot as plt
import sys
import numpy as np
from utilities import parseDAT, parseCNF


def makeDAT(datfile, files, optima, times):
    with open(datfile, 'w') as f:
        for i, file in enumerate(files):
            line = file + "   O = {0}   T = {1}\n".format(optima[i], times[i])
            f.write(line)


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print "Usage: ./categorizeDAT.py <DAT> tag"
        sys.exit(1)

    datfile = sys.argv[1]
    tag = sys.argv[2]

    files, optima, times = parseDAT(datfile)

    ratios = []

    for cnf in files:
        ratios.append(parseCNF(cnf))


    hist, edges = np.histogram(ratios, len(files))

    avgs = (edges[1:] + edges[:-1])/2

    for i in range(len(hist)):
        if hist[i] > 0:
            if i == len(hist)-1:
                indices = np.flatnonzero(np.logical_and(ratios >= edges[i], ratios <= edges[i+1]))
            else:
                indices = np.flatnonzero(np.logical_and(ratios >= edges[i], ratios < edges[i+1]))
            datfile = tag + ".{0:.3f}.dat".format(avgs[i])

            selected_files = []
            selected_optima = []
            selected_times = []

            for j in indices:
                selected_files.append(files[j])
                selected_optima.append(optima[j])
                selected_times.append(times[j])
                print("Categorized file {0} under ratio {1} in DAT: {2}".format(files[j], ratios[j], datfile))


            makeDAT(datfile, selected_files, selected_optima, selected_times)

