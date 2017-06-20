#!/usr/bin/python

import sys
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

    selected_files = {}
    selected_optima = {}
    selected_times = {}
    for i in range(len(files)):
        datfile = tag + ".{0:.3f}.dat".format(ratios[i])

        if datfile not in selected_files.keys():
            selected_files[datfile] = []
            selected_optima[datfile] = []
            selected_times[datfile] = []

        selected_files[datfile].append(files[i])
        selected_optima[datfile].append(optima[i])
        selected_times[datfile].append(times[i])

        print("Categorized file {0} under ratio {1:.3f} in DAT: {2}".format(files[i], ratios[i], datfile))


    for datfile in selected_files.keys():
        makeDAT(datfile, selected_files[datfile], selected_optima[datfile], selected_times[datfile])

