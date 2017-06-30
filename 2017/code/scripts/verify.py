#!/usr/bin/python

import sys
from subprocess import check_output

def parseSSMC(output):
    varmap = {}
    print(output)
    for line in output.split('\n'):
	if(line.startswith("v")): varmap = line.split();
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
                    break
                if var == vars[-1]:
                    sum += int(wStr)

            line = f.readline()

    return sum

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: ./analyzeSSMC.py <lut> <instance.cnf> optimum")
        sys.exit(1)

    lut = sys.argv[1]
    wcnf = sys.argv[2]
    optimal = sys.argv[3]

    args = []

    args.append("../../bin/ssmc")
    args.append(lut)
    args.append(wcnf)
    args.append(optimal)

    output = check_output(args)

    varmap = parseSSMC(output)

    sum = parseCNF(wcnf, varmap)

    print sum
