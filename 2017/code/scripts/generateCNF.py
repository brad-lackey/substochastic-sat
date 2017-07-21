#!/usr/bin/python
import sys
import numpy as np


def printToWCNF(file, lines, vars, topweight):
    # list of soft clauses' weights
    soft = [int(line.split()[0]) for line in lines if int(line.split()[0]) < topweight]
    soft_sum = sum(soft)

    if file:
        with open(file, 'w') as f:
            f.write("c Soft clauses weight: {0}\n".format(soft_sum))
            f.write("p wcnf {0} {1} {2}\n".format(vars, len(lines), topweight))
        with open(file, 'a') as f:
            for line in lines:
                f.write("{0}\n".format(line))
        print("{0} lines written to {1}.".format(len(lines) + 1, file))
    else:
        print("c Soft clauses weight: {0}".format(soft_sum))
        print("p wcnf {0} {1} {2}".format(vars, len(lines), topweight))
        for line in lines:
            print(line)


def coupleToField(lines, weight, var):
    if weight > 0:
        lines.append("{0} -{1} 0".format(weight, var))
        # print("{0} is true, cost:{1}".format(var, weight))
    else:
        lines.append("{0} {1} 0".format(-weight, var))
        # print("{0} is false, cost:{1}".format(var, -weight))


def addCoupling(lines, weight, var1, var2):
    if weight > 0:
        lines.append("{0} {1} -{2} 0".format(weight, var1, var2))
        lines.append("{0} -{1} {2} 0".format(weight, var1, var2))
        # print("{0} and {1} are same, cost: {2}".format(var1, var2, weight))
    else:
        lines.append("{0} {1} {2} 0".format(-weight, var1, var2))
        lines.append("{0} -{1} -{2} 0".format(-weight, var1, var2))
        # print("{0} and {1} are opposite, cost: {2}".format(var1, var2, -weight))


def generateCNF(clusters, Jmatrix, hvector, outfile=None):

    ws_cluster_pairs = clusters / 2

    extra_cluster = (clusters % 2 == 1)

    rows, cols, clusters_per_ws_pair = np.shape(Jmatrix)

    # number of qubits per cluster
    N = rows

    vars = ws_cluster_pairs*2*N
    if extra_cluster and ws_cluster_pairs > 1:
        vars += N

    topweight = max([2*np.max(Jmatrix), 2*np.max(hvector)])

    # ferromagnetism = 0 or 1
    fm = 0
    fm2 = 0
    lines = []
    for ws_cluster in range(ws_cluster_pairs):
        # print("---")

        for cluster in range(clusters_per_ws_pair):

            if cluster > 0 and ws_cluster > 0:
                # flip coin to determine ferromagnetism
                fm = np.random.randint(2)

                if extra_cluster and ws_cluster == ws_cluster_pairs-1:
                    fm2 = np.random.randint(2)

            for row in range(rows):
                for col in range(cols):
                    if Jmatrix[row, col, cluster]:
                        if col > row:
                            # print("2 intracluster ({0} {1})".format(N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1))
                            # couple qubits via intracluster interactions
                            addCoupling(lines, 2*Jmatrix[row,col,cluster], N*cluster + 2*N*ws_cluster + row+1, N*cluster + 2*N*ws_cluster + col+1)

                            if extra_cluster and cluster == 0 and ws_cluster_pairs > 1 and ws_cluster == ws_cluster_pairs-1:
                                # print("2 intracluster ({0} {1})".format(2*N*ws_cluster_pairs + row+1, 2*N*ws_cluster_pairs + col+1))
                                addCoupling(lines, 2*Jmatrix[row, col, cluster], 2*N*ws_cluster_pairs + row+1, 2*N*ws_cluster_pairs + col+1)

                if Jmatrix[row, row, cluster] and row+1 > N/2:
                    # print("2 intercluster (weak-strong) ({0} {1})".format(2*N*ws_cluster + row+1, 2*N*ws_cluster + N + row+1))
                    # couple qubits from cluster 1 (weak) to cluster 2 (strong) via intercluster interactions
                    addCoupling(lines, 2*Jmatrix[row,row,cluster], 2*N*ws_cluster + row+1, 2*N*ws_cluster + N + row+1)

                    # if not the first weak-strong cluster pair, then couple to previous cluster pair
                    if ws_cluster > 0:

                        if fm:
                            J = 100
                        else:
                            J = -100

                        # print("2 intercluster (strong-strong) ({0} {1})".format(2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1))
                        # couple strong clusters together
                        addCoupling(lines, 2*J, 2*N*(ws_cluster-1) + N + row+1, 2*N*ws_cluster + N + row+1)

                        # if last cluster, couple single N-qubit cluster to the rest of the pack
                        if ws_cluster_pairs > 1 and ws_cluster == ws_cluster_pairs-1:

                            if fm2:
                                J = 100
                            else:
                                J = -100

                            if extra_cluster:

                                # print("2 intercluster (strong-strong) ({0} {1})".format(2*N*(ws_cluster_pairs) + row+1, 2*N*ws_cluster + N + row+1))
                                # couple strong clusters of the last cluster to this cluster
                                addCoupling(lines, 2*J, 2*N*(ws_cluster_pairs) + row+1, 2*N*ws_cluster + N + row+1)

                                # print("2 intercluster (strong-strong) ({0} {1})".format(2*N*(ws_cluster_pairs) + row+1, N + row+1))
                                # couple strong clusters of the first cluster to this cluster
                                addCoupling(lines, 2*J, 2*N*(ws_cluster_pairs) + row+1, N + row+1)
                            elif clusters > 4:
                                # print("2 intercluster (strong-strong) ({0} {1})".format(2*N*ws_cluster + N + row+1, N + row+1))
                                # couple strong clusters of the last cluster to the first cluster
                                addCoupling(lines, 2*J, 2*N*ws_cluster + N + row+1, N + row+1)
                if cluster < 2:
                    # print("1 field-cluster ({0})".format(2*N*ws_cluster + N*cluster + row+1))
                    # couple qubits to the local field in the cluster
                    coupleToField(lines, 2*hvector[row, cluster], 2*N*ws_cluster + N*cluster + row+1)

                    # if last cluster, couple single N-qubit cluster to the field
                    if extra_cluster and cluster == 1 and ws_cluster_pairs > 1 and ws_cluster == ws_cluster_pairs-1:
                        # print("1 field-cluster ({0})".format(2*N*ws_cluster_pairs + row+1))
                        coupleToField(lines, 2*hvector[row, cluster], 2*N*ws_cluster_pairs + row+1)

    printToWCNF(outfile, lines, vars, topweight)


if __name__ == "__main__":
    if len(sys.argv) < 6 or len(sys.argv) > 7:
        print("Usage: ./generateCNF.py num_clusters qubits_per_cluster J h1(weak coupling) h2(strong coupling) [<output.wcnf>]")
        sys.exit(1)

    clusters = int(sys.argv[1])
    N = int(sys.argv[2])
    J = int(sys.argv[3])
    h1 = int(sys.argv[4])
    h2 = int(sys.argv[5])

    outfile = None
    if len(sys.argv) == 7:
        outfile = sys.argv[6]

    X = np.array([[0, 1],[1, 0]], dtype=np.int)
    Jmatrix = J*np.ones((N/2, N/2), dtype=np.int)
    Jmatrix = np.kron(X, Jmatrix)
    Jmatrix = np.dstack((Jmatrix, Jmatrix))
    Jmatrix = np.dstack((Jmatrix, np.zeros((N,N), dtype=np.int)))

    hvector = np.ones((N, 2), dtype=np.int)
    hvector[:,0] = h1
    hvector[:,1] = h2

    for i in range(N):
        Jmatrix[i, i, 0:2] = 0
        Jmatrix[i, i, 2] = J

    generateCNF(clusters, Jmatrix, hvector, outfile)

    sys.exit(0)