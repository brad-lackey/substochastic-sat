#!/usr/bin/python
import sys
import numpy as np


def generateCNF(Jmatrix, hvector, outfile=None):

    rows, cols = np.shape(Jmatrix)

    vars = rows
    cls = np.count_nonzero(Jmatrix) + rows

    topweight = max([2*np.max(Jmatrix), 2*np.max(hvector)])

    if outfile:
        with open(outfile, 'w') as f:

            f.write("p wcnf {0} {1} {2}\n".format(vars, cls, topweight))

            for row in range(rows):
                for col in range(cols):
                    if Jmatrix[row, col] and col > row:
                        f.write("{0} {1} -{2} 0\n".format(2*Jmatrix[row,col], row+1, col+1))
                        f.write("{0} -{1} {2} 0\n".format(2*Jmatrix[row,col], row+1, col+1))

                f.write("{0} -{1} 0\n".format(2*hvector[row], row+1))
    else:
        print("p wcnf {0} {1} {2}".format(vars, cls, topweight))

        for row in range(rows):
            for col in range(cols):
                if Jmatrix[row, col] and col > row:
                    print("{0} {1} -{2} 0".format(2*Jmatrix[row,col], row+1, col+1))
                    print("{0} -{1} {2} 0".format(2*Jmatrix[row,col], row+1, col+1))

            print("{0} -{1} 0".format(2*hvector[row], row+1))



if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: ./generateCNF.py N J h [<output.wcnf>]")
        sys.exit(1)

    N = int(sys.argv[1])
    J = float(sys.argv[2])
    h = float(sys.argv[3])

    outfile = None
    if len(sys.argv) == 5:
        outfile = sys.argv[4]

    Jmatrix = np.zeros((N, N))
    hvector = h*np.ones(N)

    for i in range(N-1):
        Jmatrix[i, i+1] = J
        Jmatrix[i+1, i] = J

    generateCNF(Jmatrix, hvector, outfile)

    sys.exit(0)