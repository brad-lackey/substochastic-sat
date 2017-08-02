#!/usr/bin/python
import sys
import numpy as np
from subprocess import check_output, check_call
import matplotlib.pyplot as plt
import time as TIME
from utilities import sendEmail


def parseOutput(outStr, program):
    lines = outStr.split('\n')

    time = -1
    optimum = 1e10
    if program.lstrip("./") == "ssmc":
        for line in lines:
            if line.startswith("c Walltime:"):
                time = float(line.split()[2])

            elif line.endswith("sys"):
                time = float(line.split()[0])
    elif program.lstrip("./") == "CCEHC":
        for line in lines:
            if line.startswith("c time"):
                time = float(line.split()[4])
            elif line.startswith("o"):
                optimum = int(line.split()[1])

    return time, optimum


def timeProgram(program, vars):

    # create WCNF for the given number of variables

    args = []

    args.append("./generateCNF.py")
    args.append(str(vars))
    args.append(str(8))     # qubits per cluster
    args.append(str(100))   # J = 100
    args.append(str(44))    # h1 = 44
    args.append(str(-100))  # h2 = -100
    args.append("out.wcnf")

    check_call(args)

    # solve the WCNF

    args = []

    args.append(program)
    args.append("out.wcnf")

    if program == "../../../../CCEHC":
        args.append(str(TIME.time()))
        args.append("3600")           # 1 hour of runtime

    try:
        outStr = check_output(args)
        return parseOutput(outStr, program)

    except Exception:
        print("Could not run \"{0}\"!".format(program))
        return None


def printResultsToCSV(csv_file, times, optima, N):

    with open(csv_file, 'w') as f:
        f.write("Size (Qubits)," + ",".join(map(str, N)) + "\n")
        for pgm in times.keys():
            tmp = pgm.lstrip("./") + " times,"
            tmp = tmp + ",".join(map(str, times[pgm]))
            f.write(tmp + "\n")
            tmp = pgm.lstrip("./") + " optimum,"
            tmp = tmp + ",".join(map(str, optima[pgm]))
            f.write(tmp + "\n")



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <results.csv>")
        sys.exit(1)

    csv_file = sys.argv[1]

    programs = []

    programs.append("../../../../CCEHC")

    N = [45, 64, 125]

    trials = 10

    times = {}
    optima = {}
    min_optima = {}

    for program in programs:
        for trial in range(trials):
            if program in times:
                for i, n in enumerate(N):
                    t, opt = timeProgram(program, n)
                    times[program][i] += t
                    optima[program][i].append(opt)
            else:
                for i, n in enumerate(N):
                    t, opt = timeProgram(program, n)
                    times[program] = []
                    optima[program] = []

                    times[program].append(t)
                    optima[program].append([opt])

        # compute average of times
        times[program] = [time / float(trials) for time in times[program]]

        # find the minimum optimum
        min_optima[program] = [min(opts) for opts in optima[program]]

        sendEmail("Finished analyzing {0}.".format(program.lstrip("./")))

    printResultsToCSV(csv_file, times, min_optima, N)

