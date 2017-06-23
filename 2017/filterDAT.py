#!/usr/bin/python

import sys
import numpy as np
from utilities import parseDAT, parseOUT, parseLUT
from categorizeDAT import makeDAT
from optimizeLUT import tryLUT


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print "Usage: ./filterDAT.py <LUT> <DAT> tag"
        sys.exit(1)

    lutfile = sys.argv[1]
    datfile = sys.argv[2]
    tag = sys.argv[3]

    files, optima, times = parseDAT(datfile)

    opt_dict = {}
    t_dict = {}
    for i, fname in enumerate(files):
        opt_dict[fname] = optima[i]
        t_dict[fname] = times[i]

    bins, dT, A, psize = parseLUT(lutfile)

    tryLUT("A", tag, datfile, 1, dT, A, psize)

    _files, _optima, _, _, _ = parseOUT(tag + ".out")

    hard_files = []
    hard_optima = []
    hard_times = []

    easy_files = []
    easy_optima = []
    easy_times = []

    for i, fname in enumerate(_files):

        # Filter files based on success or not

        if opt_dict[fname] != _optima[i]:
            hard_files.append(fname)
            hard_optima.append(opt_dict[fname])
            hard_times.append(t_dict[fname])
        else:
            easy_files.append(fname)
            easy_optima.append(opt_dict[fname])
            easy_times.append(t_dict[fname])

    # Create respective DAT files
    if len(easy_files) > 0:
        makeDAT(tag + ".easy.dat", easy_files, easy_optima, easy_times)
    if len(hard_files) > 0:
        makeDAT(tag + ".hard.dat", hard_files, hard_optima, hard_times)


