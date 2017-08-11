#!/usr/bin/python
import sys
import numpy as np
from subprocess32 import check_output, check_call, TimeoutExpired
import matplotlib.pyplot as plt
import time as TIME
from joblib import Parallel, delayed
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

    TIMEOUT = 3600

    if program == "../../../../CCEHC":
        args.append(str(TIME.time()))
        args.append(str(TIMEOUT))           # 1 hour of runtime

    try:
        outStr = check_output(args, timeout=1.1*TIMEOUT)
        return parseOutput(outStr, program)
    except TimeoutExpired:
        print("\"{0}\" timed out.".format(program))
        return None, None
    except Exception:
        print("Could not run \"{0}\"!".format(program))
        return None, None


def printResultsToCSV(csv_file, times, optima, avg_times, min_optima, N):

    with open(csv_file, 'w') as f:
        f.write("Size (Qubits)," + ",".join(map(str, N)) + "\n")
        for pgm in times.keys():
            tmp = pgm.lstrip("./") + " times,"
            tmp = tmp + ",".join(map(lambda x: ";".join(map(str, x)), times[pgm]))
            f.write(tmp + "\n")
            tmp = pgm.lstrip("./") + " optima,"
            tmp = tmp + ",".join(map(lambda x: ";".join(map(str, x)), optima[pgm]))
            f.write(tmp + "\n")
            tmp = pgm.lstrip("./") + " avg times,"
            tmp = tmp + ",".join(map(str, avg_times[pgm]))
            f.write(tmp + "\n")
            tmp = pgm.lstrip("./") + " min optimum,"
            tmp = tmp + ",".join(map(str, min_optima[pgm]))
            f.write(tmp + "\n")



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + " <results.csv>")
        sys.exit(1)

    csv_file = sys.argv[1]

    programs = []

    programs.append("../../../../CCEHC")

    N = [20, 25, 30, 35, 39, 45, 50, 55, 60, 64, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125]

    N_JOBS = 2

    trials = 10

    times = {}
    optima = {}
    avg_times = {}
    min_optima = {}

    for program in programs:
        full_results = Parallel(n_jobs=N_JOBS)(delayed(timeProgram)(program, n) for trial in range(trials) for n in N)
        for trial in range(trials):
            results = full_results[trial*len(N):(trial+1)*len(N)]
            for i, res in enumerate(results):
                t, opt = res
                if trial > 0:
                    times[program][i].append(t)
                    optima[program][i].append(opt)
                else:
                    if program not in times:
                        times[program] = []
                    if program not in optima:
                        optima[program] = []

                    times[program].append([t])
                    optima[program].append([opt])

        # compute average of times
        avg_times[program] = [sum(filter(None, time))/len(list(filter(None, time))) if len(list(filter(None, time))) else None for time in times[program]]

        # print the optima found
        print("Optima for {0}: {1}".format(program, optima[program]))

        # find the minimum optimum
        min_optima[program] = [min(filter(None, opts)) if len(list(filter(None, opts))) else None for opts in optima[program]]

        sendEmail("Finished analyzing {0}.".format(program.lstrip("./")))

    printResultsToCSV(csv_file, times, optima, avg_times, min_optima, N)

