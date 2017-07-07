#!/usr/bin/python
import sys
import numpy as np
from subprocess import check_output, check_call
import matplotlib.pyplot as plt

def parseOutput(outStr):
    lines = outStr.split('\n')

    for line in lines:
        if line.startswith("c Walltime:"):
            time = float(line.split()[2])

        elif line.endswith("sys"):
            time = float(line.split()[0])

    return time


def timeProgram(program, vars):

    args = []

    args.append("./generateCNF.py")
    args.append(str(vars))
    args.append(str(1))
    args.append(str(44))
    args.append(str(-100))
    args.append("out.wcnf")

    check_call(args)

    args = []

    args.append("time")
    args.append(program)
    args.append("out.wcnf")

    try:
        outStr = check_output(args)
        return parseOutput(outStr)

    except Exception:
        pass


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: ./CNFplotter.py command")
        sys.exit(1)

    program = sys.argv[1]

    N = (np.arange(10) + 1) * 100

    times = [timeProgram(program, n) for n in N]

    plt.plot(N, times)
    plt.title("Times vs. Problem-size")
    plt.xlabel("Problem size N")
    plt.ylabel("Time in seconds to completion")
    plt.show()

    print times, N