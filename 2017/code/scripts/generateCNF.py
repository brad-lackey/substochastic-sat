#!/usr/bin/python
import sys
import numpy as np


def generateCNF(Jmatrix, hvector, outfile=None):

    # number of weak-strong cluster pairs
    ws_cluster_pairs = 19

    rows, cols, clusters = np.shape(Jmatrix)

    # number of qubits per cluster
    N = rows

    vars = ws_cluster_pairs*2*N

    if N == 1:
        cls = 4
    else:
        ws_clauses = (np.count_nonzero(Jmatrix[:,:,0:2]) + 2*rows + rows)
        cls = ws_cluster_pairs * ws_clauses

        if ws_cluster_pairs > 1:
            cls += N * (ws_cluster_pairs-1)
            cls += 2*N


    topweight = max([2*np.max(Jmatrix), 2*np.max(hvector)])

    if outfile:
        with open(outfile, 'w') as f:
            f.write("p wcnf {0} {1} {2}\n".format(vars, cls, topweight))
    else:
        print("p wcnf {0} {1} {2}".format(vars, cls, topweight))

    lines = 1
    for ws_cluster in range(ws_cluster_pairs):
        # print("---")
        if outfile:
            with open(outfile, 'a') as f:

                for row in range(rows):
                    for cluster in range(clusters):
                        for col in range(cols):
                            if Jmatrix[row, col, cluster]:
                                if col > row:
                                    # print("2 intracluster ({0} {1})".format(N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1))
                                    lines += 2
                                    # couple qubits via intracluster interactions
                                    f.write("{0} {1} -{2} 0\n".format(2*Jmatrix[row,col,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1))
                                    f.write("{0} -{1} {2} 0\n".format(2*Jmatrix[row,col,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1))
                        if Jmatrix[row, row, cluster] and row+1 > N/2:
                            # print("2 intercluster (weak-strong) ({0} {1})".format(N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + N + row+1))
                            lines += 2
                            # couple qubits from cluster 1 (weak) to cluster 2 (strong) via intercluster interactions
                            f.write("{0} {1} -{2} 0\n".format(2*Jmatrix[row,row,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + N + row+1))
                            f.write("{0} -{1} {2} 0\n".format(2*Jmatrix[row,row,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + N + row+1))

                            # if not the first weak-strong cluster pair, then couple to previous cluster pair
                            if ws_cluster > 0:
                                # flip coin to determine ferromagnetism
                                coin = np.random.randint(2)

                                if coin == 0:
                                    J = 1
                                else:
                                    J = -1

                                # print("2 intercluster (strong-strong) ({0} {1})".format(2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1))
                                lines += 2
                                # couple strong clusters together
                                f.write("{0} {1} -{2} 0\n".format(2*J, 2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1))
                                f.write("{0} -{1} {2} 0\n".format(2*J, 2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1))

                                # if last cluster, couple single N-qubit cluster to the rest of the pack
                                if ws_cluster_pairs > 1 and ws_cluster == ws_cluster_pairs-1:

                                    # flip coin to determine ferromagnetism
                                    coin = np.random.randint(2)

                                    if coin == 0:
                                        J = 1
                                    else:
                                        J = -1

                                    # print("2 intercluster (strong-strong) ({0} {1})".format( 2*N*(ws_cluster_pairs) + N + row+1, 2*N*ws_cluster + N + row+1))
                                    # couple strong clusters of the last cluster to this cluster
                                    f.write("{0} {1} -{2}\n".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, 2*N*ws_cluster + N + row+1))
                                    f.write("{0} -{1} {2}\n".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, 2*N*ws_cluster + N + row+1))

                                    # print("2 intercluster (strong-strong) ({0} {1})".format(2*N*(ws_cluster_pairs) + N + row+1, N + row+1))
                                    # couple strong clusters of the first cluster to this cluster
                                    f.write("{0} {1} -{2}\n".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, N + row+1))
                                    f.write("{0} -{1} {2}\n".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, N + row+1))

                                    lines += 4


                        if cluster < 2:
                            # print("1 field-cluster ({0})".format(2*N*ws_cluster + N*cluster + row+1))
                            lines += 1
                            # couple qubits to the local field in the cluster
                            f.write("{0} -{1} 0\n".format(2*hvector[row, cluster], 2*N*ws_cluster + N*cluster + row+1))

        else:

            for row in range(rows):
                for cluster in range(clusters):
                    for col in range(cols):
                        if Jmatrix[row, col, cluster]:
                            if col > row:
                                # couple qubits via intracluster interactions
                                print("{0} {1} -{2} 0".format(2*Jmatrix[row,col,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1))
                                print("{0} -{1} {2} 0".format(2*Jmatrix[row,col,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1))
                    if Jmatrix[row, row, cluster] and row+1 >= N/2:
                        # couple qubits via intercluster interactions
                        print("{0} {1} -{2} 0".format(2*Jmatrix[row,row,cluster], 2*N*ws_cluster + row+1, 2*N*ws_cluster + N + row+1))
                        print("{0} -{1} {2} 0".format(2*Jmatrix[row,row,cluster], 2*N*ws_cluster + row+1, 2*N*ws_cluster + N + row+1))


                        # if not the first weak-strong cluster pair, then couple to previous cluster pair
                        if ws_cluster > 0:
                            # flip coin to determine ferromagnetism
                            coin = np.random.randint(2)

                            if coin == 0:
                                J = 1
                            else:
                                J = -1

                            # couple strong clusters together
                            print("{0} {1} -{2} 0".format(2*J, 2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1))
                            print("{0} -{1} {2} 0".format(2*J, 2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1))

                            # if last cluster, couple single N-qubit cluster to the rest of the pack
                            if ws_cluster_pairs > 1 and ws_cluster == ws_cluster_pairs-1:

                                # flip coin to determine ferromagnetism
                                coin = np.random.randint(2)

                                if coin == 0:
                                    J = 1
                                else:
                                    J = -1

                                # couple strong clusters of the last cluster to this cluster
                                print("{0} {1} -{2}".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, 2*N*ws_cluster + N + row+1))
                                print("{0} -{1} {2}".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, 2*N*ws_cluster + N + row+1))

                                # couple strong clusters of the first cluster to this cluster
                                print("{0} {1} -{2}".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, N + row+1))
                                print("{0} -{1} {2}".format(2*J, 2*N*(ws_cluster_pairs) + N + row+1, N + row+1))


                    if cluster < 2:
                        # couple qubits to the local field in the cluster
                        print("{0} -{1} 0".format(2*hvector[row, cluster], 2*N*ws_cluster + N*cluster + row+1))

    if outfile:
        print("{0} lines written to {1}".format(lines, outfile))


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