#!/usr/bin/python

import sys
from subprocess import check_output

def parseSSMC(output):
    varmap = {}
    print(output)
    for line in output.split('\n'):
        if(line.startswith("v")):
            varmap = line.split()

        return varmap


def parseCNF(cnf, varmap):
    sum = 0

    with open(cnf, 'r') as f:
        line = f.readline()
        while(len(line) > 0):
            if line.startswith("p") or line.startswith("c"):
                line = f.readline()
                continue

            vars = line.split()

            wStr = vars[0]
            for var in vars[1:-1]:
                if var in varmap:
                    # pass
                    break
                if var == vars[-2]:
                    sum += int(wStr)

            line = f.readline()

    return sum

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage: ./verify.py <instance.cnf> x1 x2 .. xN")
        sys.exit(1)

    wcnf = sys.argv[1]
    varmap = sys.argv[2:]

    sum = parseCNF(wcnf, varmap)

    print sum
