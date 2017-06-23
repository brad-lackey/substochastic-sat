#!/usr/bin/python

import sys
from utilities import parseDAT, parseOUT, parseLUT, parseCNF
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
    i = 0

    while i < len(files):
        cnf = files[i]
        variables, clauses = parseCNF(cnf)

        if variables >= 500:
            del files[i]
            del optima[i]
            del times[i]
        else:
            opt_dict[cnf] = optima[i]
            t_dict[cnf] = times[i]
            i += 1

    datfile = tag + ".dat"

    if len(files) > 0:
        makeDAT(datfile, files, optima, times)

    bins, dT, A, psize = parseLUT(lutfile)

    tryLUT("A", tag, datfile, 1, dT, A, psize)

    _files, _optima, _, _, _ = parseOUT(tag + ".out")

    hard_files = []
    hard_optima = []
    hard_times = []

    easy_files = []
    easy_optima = []
    easy_times = []

    for i, cnf in enumerate(_files):

        # Filter files based on success or not

        if opt_dict[cnf] != _optima[i]:
            hard_files.append(cnf)
            hard_optima.append(opt_dict[cnf])
            hard_times.append(t_dict[cnf])
        else:
            easy_files.append(cnf)
            easy_optima.append(opt_dict[cnf])
            easy_times.append(t_dict[cnf])

    # Create respective DAT files
    if len(easy_files) > 0:
        makeDAT(tag + ".easy.dat", easy_files, easy_optima, easy_times)
    if len(hard_files) > 0:
        makeDAT(tag + ".hard.dat", hard_files, hard_optima, hard_times)


