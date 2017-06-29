#!/usr/bin/python

import sys
from subprocess import check_output

def parseSSMC(output):
    varmap = {}
    for line in output.split('\n'):
        if line.startswith("v"):
            for v in line.split():
                varmap[v] = not v.startswith('-')

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

            if len(vars) > 3:
                pass

            for var in vars[1:-1]:
                if(var not in varmap):
                    if var == vars[-2]:
                        sum += int(wStr)
                else:
                    break

            line = f.readline()

    return sum

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: ./analyzeSSMC.py <lut> <instance.cnf> optimum")

    lut = sys.argv[1]
    wcnf = sys.argv[2]
    optimal = sys.argv[3]

    args = []

    args.append("./ssmc")
    args.append(lut)
    args.append(wcnf)
    args.append(optimal)

    output = check_output(args)

    varmap = parseSSMC(output)

    sum = parseCNF(wcnf, varmap)

    print sum