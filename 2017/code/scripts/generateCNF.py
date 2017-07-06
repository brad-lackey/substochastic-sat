#!/usr/bin/python
import sys
import numpy as np


def generateCNF(Jmatrix, hvector, outfile=None):

    rows, cols, clusters = np.shape(Jmatrix)

    N = rows

    vars = 2*rows
    cls = np.count_nonzero(Jmatrix[:,:,0:2]) + 2*N + N + 2

    topweight = max([2*np.max(Jmatrix), 2*np.max(hvector)])

    if outfile:
        with open(outfile, 'w') as f:
            f.write("p wcnf {0} {1} {2}\n".format(vars, cls, topweight))

            for row in range(rows):
                for cluster in range(clusters):
                    for col in range(cols):
                        if Jmatrix[row, col, cluster]:
                            if col > row:
                                f.write("{0} {1} -{2} 0\n".format(2*Jmatrix[row,col,cluster], N*cluster + row+1, N*cluster + col+1))
                                f.write("{0} -{1} {2} 0\n".format(2*Jmatrix[row,col,cluster], N*cluster + row+1, N*cluster + col+1))
                    if Jmatrix[row, row, cluster] and row+1 >= N/2:
                        f.write("{0} {1} -{2} 0\n".format(2*Jmatrix[row,row,cluster], row+1, N + row+1))
                        f.write("{0} -{1} {2} 0\n".format(2*Jmatrix[row,row,cluster], row+1, N + row+1))


                    if cluster < 2:
                        f.write("{0} -{1} 0\n".format(2*hvector[row, cluster], N*cluster + row+1))

        print("WCNF written to {0}".format(outfile))
    else:
        print("p wcnf {0} {1} {2}".format(vars, cls, topweight))

        for row in range(rows):
            for cluster in range(clusters):
                for col in range(cols):
                    if Jmatrix[row, col, cluster]:
                        if col > row:
                            print("{0} {1} -{2} 0".format(2*Jmatrix[row,col,cluster], N*cluster + row+1, N*cluster + col+1))
                            print("{0} -{1} {2} 0".format(2*Jmatrix[row,col,cluster], N*cluster + row+1, N*cluster + col+1))
                if Jmatrix[row, row, cluster] and row+1 >= N/2:
                    print("{0} {1} -{2} 0".format(2*Jmatrix[row,row,cluster], row+1, N + row+1))
                    print("{0} -{1} {2} 0".format(2*Jmatrix[row,row,cluster], row+1, N + row+1))


                if cluster < 2:
                    print("{0} -{1} 0".format(2*hvector[row, cluster], N*cluster + row+1))



if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Usage: ./generateCNF.py N J h1 h2 [<output.wcnf>]")
        sys.exit(1)

    N = int(sys.argv[1])
    J = int(sys.argv[2])
    h1 = int(sys.argv[3])
    h2 = int(sys.argv[4])

    outfile = None
    if len(sys.argv) == 6:
        outfile = sys.argv[5]

    Jmatrix = J*np.ones((N, N, 3), dtype=np.int)
    Jmatrix[:,:,2] = 0
    hvector = np.ones((N, 2), dtype=np.int)
    hvector[:,0] = h1
    hvector[:,1] = h2

    for i in range(N):
        Jmatrix[i, i, 0:2] = 0
        Jmatrix[i, i, 2] = J


    generateCNF(Jmatrix, hvector, outfile)

    sys.exit(0)