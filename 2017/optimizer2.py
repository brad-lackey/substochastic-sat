#!/usr/bin/python

from optimizeLUT import parseLUT, sendEmail, getABounds, plotLUT, plotPsize
from scipy import stats
from subprocess32 import check_call, TimeoutExpired
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fminbound
import datetime
import sys
from createLUT import makeLUT
from histAnalysis import parseOUT

BOUND_CAP = 0.1  # cap on the bounds
BOUND_MULTIPLIER = 1.1  # fraction over which the bound can extend
UPDATE_PENALTY = 10000000  # penalty to give scripts which timeout
N_ITERS_CAP = 5  # max number of optimization iterations
RECURSION_LIMIT = 5  # max levels optimizer can branch LUT
THRESHOLD = 0.1  # min threshold before accepting new minimum

def main():
    args = sys.argv
    email = False
    verbose = False
    plotenabled = False

    if '-m' in args:
        email = True
        args.remove('-m')
    if '-v' in args:
        verbose = True
        args.remove('-v')
    if '-p' in args:
        plotenabled = True
        args.remove('-p')
    if len(args) == 7 or len(args) == 9:
        global var
        var = args[1]
        global xpmt
        xpmt = int(args[2])

        if xpmt < 0 or xpmt > 2:
            print("Usage: ./optimizer2 dT|A|psize|both <experiment_type (0:fwd/bwd, 1:2 rnd, 2:2 no-cons rnd)> [-v] [-m] [-p] <initialLUT> <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
            return 1

        lutfile = args[3]
        datfile = args[4]
        trials = args[5]
        tag = args[6]
        weight = None
        runtime = None

        if len(args) == 9:
            weight = args[7]
            runtime = args[8]

    else:
        print("Usage: ./optimizer2 dT|A|psize|both <experiment_type (0:fwd/bwd, 1:2 rnd, 2:2 no-cons rnd)> [-v] [-m] [-p] <initialLUT> <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
        return 1

    optimizeLUT(var, lutfile, datfile, trials, tag, weight, runtime, recursion_level=0, email=email, verbose=verbose, plotenabled=plotenabled)

    return 0


"""Returns the negative t-statistic from a t-test from this set of conf files using given LUT to the previous best LUT"""
def tryLUT(tag, filename, trials, dT, A, psize, weight=None, runtime=None, plotenabled=False, verbose=False):
    if len(dT) != len(A) or len(psize) != len(dT) or len(psize) != len(A):
        raise Exception("Vectors dT, A and psize are not the same length!")

    bins = len(dT)
    lut = tag + ".LUT.txt"

    psize = map(round, psize)
    psize = map(int, psize)

    makeLUT(lut, bins, dT, A, psize)

    args = []
    args.append('./testrun.pl')  # the program to run
    args.append('./ssmc')
    args.append(lut)
    args.append(filename)
    args.append(str(trials))
    args.append(tag)

    if weight and runtime:
        args.append(weight)
        args.append(runtime)

    # returns 0 if successful otherwise throws error
    try:
        check_call(args)
    except TimeoutExpired:
        return UPDATE_PENALTY

    outfile = tag + ".out"
    _, times, loops, updates = parseOUT(outfile)

    global best_updates
    if best_updates is None:
        best_updates = updates
        return 0.0
    else:
        tstat, p = stats.ttest_ind(updates, best_updates)

    if plotenabled:
        global var
        if var == 'psize':
            plotPsize(dT, psize)
        else:
            plotLUT(dT, A)

    if verbose:
        print("Tried dT=" + str(dT) + ", A=" + str(A) + ", Psize=" + str(psize) + "  with a t-stat=" + str(tstat) + ", p={0}".format(p))

    if p < THRESHOLD:
        return tstat
    else:
        return 0.0


def optimizeLUT(var, lutfile, datfile, trials, tag, weight, runtime, recursion_level=0, email=False, verbose=False, plotenabled=False, start=datetime.datetime.now()):
    if recursion_level == 0:
        if plotenabled:
            # Turn on interactive plotting
            plt.ion()
        if verbose:
            print("########## STARTING OPTIMIZATION - " + datetime.datetime.now().strftime(
                "%a %d/%m/%y %H:%M:%S") + " ##########")

    # Load initial conditions for dT and A
    bins, dT, A, psize = parseLUT(lutfile)

    f = getMinimizer(var)

    def printf(fval):
        if verbose:
            print("#################### Found new minimum: " + str(fval) + " ####################")
        if email:
            msg = "Found new minimum: " + str(fval)
            sendEmail(msg)

    # set initial minimum to initial LUT performance
    global best_updates
    best_updates = None
    fmin = tryLUT(tag, datfile, trials, dT, A, psize, weight, runtime, plotenabled, verbose)

    if recursion_level >= RECURSION_LIMIT:
        return fmin, dT, A, psize

    if var == 'dT':
        varmin = dT.copy()
    elif var == 'A':
        varmin = A.copy()
    elif var == "psize":
        varmin = psize.copy()

    if xpmt == 0:
        indices = np.concatenate((np.arange(bins), np.arange(bins-1)[::-1]))  # [0 1 2 .. bins-1 .. 2 1 0]
    else:
        indices = np.concatenate((np.arange(bins), np.arange(bins)))

    if var != "both":

        newLUT = False

        # i -> iteration
        for i in range(N_ITERS_CAP):
            fval = 0

            if xpmt != 0:
                np.random.shuffle(indices)  # shuffles the indices array for random choice of index to optimize

            minFound = False

            for row in indices:
                if var == "dT":

                    if xpmt == 1:
                        # skip for the last edge
                        if row == bins-1:
                            continue

                        varvector, othervector = dT, A

                        edges = np.insert(np.cumsum(dT), 0, 0)

                        lbound, ubound = edges[row], edges[row+2]

                        x0, fval, ierr, numfunc = fminbound(f, lbound, ubound, args=(
                            row, tag, datfile, trials, varvector, othervector, psize, weight, runtime),
                                                            full_outplotPsizeput=True, xtol=0.01)

                        edges[row+1] = x0

                        dT = np.diff(edges)
                        varvector = dT
                    else:
                        varvector, othervector = dT, A

                        lbound, ubound = 0.1, 100.0

                        x0, fval, ierr, numfunc = fminbound(f, lbound, ubound, args=(
                            row, np.delete(varvector, row), tag, datfile, trials, othervector, psize, weight,
                            runtime, verbose, plotenabled), full_output=True, xtol=0.01)
                        varvector[row] = x0

                elif var == "A":
                    varvector, othervector = A, dT

                    lbound, ubound = 0.1, 1.0

                    x0, fval, ierr, numfunc = fminbound(f, lbound, ubound, args=(
                        row, np.delete(varvector, row), tag, datfile, trials, othervector, psize, weight,
                        runtime, verbose, plotenabled), full_output=True, xtol=0.01)
                    varvector[row] = x0

                elif var == "psize":
                    varvector, othervector = psize, dT

                    lbound, ubound = 0.5 * psize[row], 2 * psize[row]

                    x0, fval, ierr, numfunc = fminbound(f, lbound, ubound, args=(
                        row, np.delete(varvector, row), tag, datfile, trials, othervector, A, weight,
                        runtime, verbose, plotenabled), full_output=True, xtol=1)
                    varvector[row] = int(round(x0))

                if fval < 0:
                    fmin = fval

                    _, _, _, best_updates = parseOUT(tag + ".out")

                    varmin = varvector.copy()

                    printf(fmin)

                    # Store the best var
                    lut = tag + ".OPTIMAL." + var + ".txt"
                    if var == "dT":
                        makeLUT(lut, bins, varmin, A, psize)
                    elif var == "psize":
                        makeLUT(lut, bins, dT, A, varmin)
                    else:
                        makeLUT(lut, bins, dT, varmin, psize)

                    if plotenabled:
                        plt.savefig(tag + ".OPTIMAL." + var + ".png")

                    newLUT = minFound = True

                if verbose:
                    print(
                        "---------- Found {0}[{1}]={2}".format(var, row, varvector[row]) + " at updates " + str(fval) + " after " + str(
                            numfunc) + " tries, {0}/{1} iterations ----------".format(i + 1, N_ITERS_CAP))

            if email:
                msg = "Progress: {0}/{1} iterations complete.".format(i + 1, N_ITERS_CAP) + "\n"
                msg += "Level: {0}\n".format(recursion_level)
                # msg += "Time spent so far: " + str(datetime.datetime.now() - start) + '\n'
                sendEmail(msg)

            if not minFound:
                print("No changes detected. Breaking out of loop...")

                msg = "Optimization converged at {0}, with {1}={2}".format(fmin, var, varmin)
                if email:
                    sendEmail(msg)
                if verbose:
                    print(msg)

                if not newLUT:
                    msg = "No improvements after {0} levels. Finishing up...".format(recursion_level)
                    if verbose:
                        print(msg)
                    if email:
                        sendEmail(msg)
                    return fmin, dT, A, psize
                break

    else:
        fmin1 = fmin2 = fmin

        changedLUT = False

        lut = tag + ".LUT.txt"
        makeLUT(lut, bins, dT, A, psize)

        while True:

            fmin1, new_dT, new_A, new_psize = optimizeLUT('A', lut, datfile, trials, tag, weight, runtime,
                                               recursion_level=recursion_level, email=email, verbose=verbose,
                                               plotenabled=plotenabled, start=start)

            if fmin1 < fmin2 and fmin1 < fmin:
                makeLUT(lut, bins, new_dT, new_A, new_psize)
                changedLUT = True

            fmin2, new_dT, new_A, new_psize = optimizeLUT('dT', lut, datfile, trials, tag, weight, runtime,
                                               recursion_level=recursion_level, email=email, verbose=verbose,
                                               plotenabled=plotenabled, start=start)

            if fmin2 < fmin1 and fmin2 < fmin:
                makeLUT(lut, bins, new_dT, new_A, new_psize)
                changedLUT = True

            if fmin1 > fmin and fmin2 > fmin:
                if changedLUT:
                    # if it cannot improve it past the fmin, save the best schedule and break the loop
                    if verbose:
                        print("Cannot improve past fmin. Breaking out of loop...")
                    break
                else:
                    # if it hasn't improved the given schedule at all, return and break out of the recursion...our job
                    #  is done.
                    if verbose:
                        print("No improvements detected. Returning fmin = {0}".format(fmin))
                    return fmin, dT, A, psize
            else:
                fmin = min([fmin1, fmin2])

    if var == "both":
        lut = tag + ".LUT.txt"
        fmin, dT, A = branchLUT(lut, tag, datfile, trials, weight, runtime, recursion_level, email, plotenabled, verbose, start)
    elif var == 'A':
        A = varmin.copy()
    elif var == 'psize':
        psize = varmin.copy()
    else:
        dT = varmin.copy()

    if recursion_level == 0:
        if verbose:
            # Print the best updates
            print("Best # updates: " + str(fmin))
        if plotenabled:
            if var == "psize":
                plotPsize(dT, psize)
            else:
                plotLUT(dT, A)
        if email:
            if var == "both":
                msg = "Optimization finished after " + str(datetime.datetime.now() - start) + \
                      ", at " + datetime.datetime.now().strftime("%a %d/%m/%y %H:%M:%S") + \
                      "!\nOptimal dT: {0}\nOptimal A: {1}\nOptimum # updates: {2}\n".format(dT, A, fmin)
            else:
                msg = "Optimization finished after " + str(datetime.datetime.now() - start) + \
                      ", at " + datetime.datetime.now().strftime("%a %d/%m/%y %H:%M:%S") + \
                      "!\nOptimal " + var + ": " + str(varmin) + "\nOptimum # updates: " + str(fmin) + "\n"
            sendEmail(msg)

    return fmin, dT, A, psize


def getMinimizer(var):
    # Minimize var
    if var == 'dT':
        if xpmt == 1:
            def f(edge, edgeI, tag, filename, trials, dT, A, psize, weight, runtime, p=False, v=False):
                edges = np.insert(np.cumsum(dT), 0, 0)

                edges[edgeI + 1] = edge

                dT = np.diff(edges)

                return tryLUT(tag, filename, trials, dT, A, psize, weight, runtime, p, v)
        else:
            f = lambda x1, i, x2, a1, a2, a3, a4, psize, a5, a6, v, p: tryLUT(a1, a2, a3, np.insert(x2, i, x1),
                                                                              a4, psize, a5,
                                                                              a6, verbose=v,
                                                                              plotenabled=p)  # rearranging the arguments for dT
    elif var == 'A':
        f = lambda x1, i, x2, a1, a2, a3, a4, psize, a5, a6, v, p: tryLUT(a1, a2, a3, a4, np.insert(x2, i, x1), psize, a5,
                                                                          a6, verbose=v,
                                                                          plotenabled=p)  # rearranging the arguments for A
    elif var == 'both':
        return None
    elif var == 'psize':
        f = lambda x1, i, x2, a1, a2, a3, a4, a5, a6, a7, v, p: tryLUT(a1, a2, a3, a4, a5, np.insert(x2, i, x1), a6, a7,
                                                                       verbose=v, plotenabled=p)
    else:
        raise Exception("Invalid variable argument! Must be \"dT\", \"A\" or \"both\"")
    return f


def branchLUT(lut, tag, datfile, trials, weight, runtime, recursion_level, email, plotenabled, verbose, start):
    # Get the best var
    bins, dT, A, psizes = parseLUT(lut)
    # split dT's largest bin in half
    maxBin = dT.max()
    maxI = dT.argmax()
    dT[maxI] = maxBin / 2.0
    dT = np.insert(dT, maxI, maxBin / 2.0)
    A = np.insert(A, maxI, A[maxI])
    psizes = np.insert(psizes, maxI, psizes[maxI])

    lut = lut.rstrip(".txt") + "." + str(recursion_level) + ".txt"

    if verbose:
        print("###########################################################")
        print("###########################################################")

    makeLUT(lut, bins + 1, dT, A, psizes)

    if verbose:
        print("Recursing down to level {0}...".format(recursion_level + 1))
        print("dT={0}\nA={1}".format(dT, A))
    if email:
        msg = "Recursing down to level {0}...".format(recursion_level+1)
        sendEmail(msg)

    return optimizeLUT("both", lut, datfile, trials, tag, weight, runtime, recursion_level=recursion_level + 1,
                       email=email, verbose=verbose, plotenabled=plotenabled, start=start)


if __name__ == "__main__":
    sys.exit(main())
