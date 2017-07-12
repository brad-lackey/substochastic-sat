#!/usr/bin/python
import sys
import numpy as np


def couple(lines, var1, var2, equal):
    if equal:
        lines.append("2 {0} -{1} 0".format(var1, var2))
        lines.append("2 -{0} {1} 0".format(var1, var2))
    else:
        lines.append("2 {0} {1} 0".format(var1, var2))
        lines.append("2 -{0} -{1} 0".format(var1, var2))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: ./bitstringToWCNF.py <bitstring>")
        sys.exit(1)

    # exclude the filename
    bitstring = sys.argv[1:]
    bits = [not bitStr.startswith('-') for bitStr in bitstring]

    # bit indices start from 1
    bitIndices = [i+1 for i, val in enumerate(bitstring)]

    topweight = 100

    lines = []
    # convert to clauses
    for i, v1 in enumerate(bits):
        # print("Var 1: {0}".format(bitIndices[i]))
        for j in (bitIndices[i+1:]):
            # print("Var 2: {0}".format(bitIndices[j-1]))
            if bits[i] == bits[j-1]:
                # bits are equal
                couple(lines, i+1, j, True)
            else:
                # bits are opposite
                couple(lines, i+1, j, False)
        # add field couplings
        lines.append("{0} {1} 0".format(topweight, bitstring[i]))

    print("c Answer: {0}, Cost: 0".format(" ".join(bitstring)))
    # print wcnf file to output
    print("p wcnf {0} {1} {2}".format(len(bits), len(lines), topweight))
    for line in lines:
        print(line)