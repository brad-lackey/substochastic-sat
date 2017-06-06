#!/usr/bin/python
import sys
from createLUT import makeLUT
from subprocess32 import check_call, TimeoutExpired
import numpy as np
from scipy.optimize import fminbound
import matplotlib.pyplot as plt
import smtplib
import datetime

BOUND_CAP = 0.1  # cap on the bounds
BOUND_MULTIPLIER = 1.1  # fraction over which the bound can extend
UPDATE_PENALTY = 100000  # penalty to give scripts which timeout
N_ITERS_CAP = 1  # max number of optimization iterations
RECURSION_LIMIT = 5  # max levels optimizer can branch LUT

"""Returns the dT and A vectors from a LUT file as a tuple"""
def parseLUT(lutfile):

    with open(lutfile, 'r') as f:
        line = f.readline()
        bins = int(line)
        dT = np.zeros(bins)
        A = np.zeros(bins)
        row = 0
        while len(line) > 0:
            line = f.readline()
            if len(line) > 0:
                dT[row], A[row] = line.rstrip('\n').split('\t')
                row += 1
        if row != bins:
            raise Exception("Invalid LUT file format!")

    return bins, dT, A

"""Returns the percentage of hits and avg runtime as a tuple"""
def parseTXT(txtfile):
    last = ''
    with open(txtfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            last = line
            line = f.readline()

    # Parse last line
    _, hitStr, _, tStr, lStr, _, uStr, _ = last.split()

    fraction = map(float, hitStr[0:hitStr.rfind('(')].split('/'))
    hit = fraction[0]/fraction[1]
    loops = float(lStr)
    t = float(tStr.rstrip("s"))
    updates = float(uStr)

    return hit, updates

"""Returns the avg updates of a set of conf files using given LUT"""
def tryLUT(tag, filename, trials, dT, A, weight=None, runtime=None, plotenabled=False, verbose=False):
    if len(dT) != len(A):
        raise Exception("Vectors dT and A are not the same length!")

    bins = len(dT)
    lut = tag + ".LUT.txt"

    makeLUT(lut, bins, dT, A)

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

    txtfile = tag + ".txt"
    hits, updates = parseTXT(txtfile)

    if hits < 1:
        updates = UPDATE_PENALTY

    if plotenabled:
        # Plot A vs t
        t = np.cumsum(dT)
        t = t - np.ediff1d(t, to_begin=t[0])/2.0  # staggers the time so that it falls in between the bins
        plt.hold(False)
        plt.plot(t, A)
        plt.ylabel("A-Values")
        plt.xlabel("Time")
        plt.title("A vs. T")
        ax = plt.gca()
        ax.relim()
        ax.autoscale_view()
        plt.draw()

    if verbose:
        print("Tried dT=" + str(dT) + ", A=" + str(A) + " with updates=" + str(updates))

    return updates


def sendEmail(msg):
    server = smtplib.SMTP("smtp.gmail.com", 587)
    server.starttls()
    server.login("email.notifier.bryanluu@gmail.com", "7788382652")
    server.sendmail("email.notifier.bryanluu@gmail.com", "bryanluu30794@gmail.com", '\n'+msg)
    server.quit()


def getABounds(bins, row, varvector):
    if row == 0 or row == bins - 1:
        if row == 0:
            diff = BOUND_MULTIPLIER * abs(varvector[row + 1] - varvector[row])
        else:
            diff = BOUND_MULTIPLIER * abs(varvector[row - 1] - varvector[row])

        if diff > BOUND_CAP:
            diff = BOUND_CAP

        lbound = varvector[row] - diff
        ubound = varvector[row] + diff
    else:
        ldiff = BOUND_MULTIPLIER * abs(varvector[row] - varvector[row - 1])
        rdiff = BOUND_MULTIPLIER * abs(varvector[row + 1] - varvector[row])

        if ldiff > BOUND_CAP:
            ldiff = BOUND_CAP
        if rdiff > BOUND_CAP:
            rdiff = BOUND_CAP

        if varvector[row - 1] <= varvector[row + 1]:
            lbound = varvector[row] - ldiff
            ubound = varvector[row] + rdiff
        else:
            ubound = varvector[row] + ldiff
            lbound = varvector[row] - rdiff
    lbound = max(lbound, 0)
    ubound = min(ubound, 1)
    return lbound, ubound


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
        var = args[1]
        global xpmt
        xpmt = int(args[2])

        if xpmt < 0 or xpmt > 2:
            print("Usage: ./optimizeLUT dT|A|both <experiment_type (0:fwd/bwd, 1:2 rnd, 2:2 no-cons rnd)> [-v] [-m] [-p] <initialLUT> <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
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
        print("Usage: ./optimizeLUT dT|A|both <experiment_type (0:fwd/bwd, 1:2 rnd, 2:2 no-cons rnd)> [-v] [-m] [-p] <initialLUT> <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
        return 1

    optimizeLUT(var, lutfile, datfile, trials, tag, weight, runtime, recursion_level=0, email=email, verbose=verbose, plotenabled=plotenabled)


    return 0


def optimizeLUT(var, lutfile, datfile, trials, tag, weight, runtime, recursion_level=0, email=False, verbose=False, plotenabled=False):
    if recursion_level == 0:
        if plotenabled:
            # Turn on interactive plotting
            plt.ion()
        if verbose:
            start = datetime.datetime.now()
            print("########## STARTING OPTIMIZATION - " + datetime.datetime.now().strftime(
                "%a %d/%m/%y %H:%M:%S") + " ##########")

    # Load initial conditions for dT and A
    bins, dT, A = parseLUT(lutfile)

    f = getMinimizer(var)

    def printf(fval):
        if verbose:
            print("#################### Found new minimum: " + str(fval) + " ####################")
        if email:
            msg = "Found new minimum: " + str(fval)
            sendEmail(msg)

    # set initial minimum to initial LUT performance
    fmin = tryLUT(tag, datfile, trials, dT, A, weight, runtime, plotenabled, verbose)

    if recursion_level >= RECURSION_LIMIT:
        return fmin, dT, A

    if var == 'dT':
        varmin = dT.copy()
    elif var == 'A':
        varmin = A.copy()

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
                            row, tag, datfile, trials, varvector, othervector, weight, runtime),
                                                            full_output=True, xtol=0.01)

                        edges[row+1] = x0

                        dT = np.diff(edges)
                        varvector = dT
                    else:
                        varvector, othervector = dT, A

                        lbound, ubound = 0.1, 2.0

                        x0, fval, ierr, numfunc = fminbound(f, lbound, ubound, args=(
                            row, np.delete(varvector, row), tag, datfile, trials, othervector, weight, runtime, verbose, plotenabled),
                                                            full_output=True, xtol=0.01)
                        varvector[row] = x0

                elif var == "A":
                    varvector, othervector = A, dT

                    lbound, ubound = getABounds(bins, row, varvector)

                    x0, fval, ierr, numfunc = fminbound(f, lbound, ubound, args=(
                    row, np.delete(varvector, row), tag, datfile, trials, othervector, weight, runtime, verbose, plotenabled),
                                                        full_output=True, xtol=0.01)
                    varvector[row] = x0

                if fval < fmin:
                    fmin = fval

                    varmin = varvector.copy()

                    printf(fmin)

                    # Store the best var
                    lut = tag + ".OPTIMAL." + var + ".txt"
                    if var == "dT":
                        makeLUT(lut, bins, varmin, A)
                    else:
                        makeLUT(lut, bins, dT, varmin)

                    if plotenabled:
                        plt.savefig(tag + ".OPTIMAL." + var + ".png")

                    newLUT = minFound = True

                if verbose:
                    print(
                    "---------- Found {0}[{1}]={2}".format(var, row, x0) + " at updates " + str(fval) + " after " + str(
                        numfunc) + " tries, {0}/{1} iterations ----------".format(i + 1, N_ITERS_CAP))

            if email:
                msg = "Progress: {0}/{1} iterations complete.".format(i + 1, N_ITERS_CAP) + "\n"
                msg += "Time spent so far: " + str(datetime.datetime.now() - start) + '\n'
                sendEmail(msg)

            if not minFound:
                print("No changes detected. Breaking out of loop...")
                msg = "Optimization converged at {0}, with {1}={2}".format(fmin, var, varmin)
                sendEmail(msg)

                if not newLUT:
                    return fmin, dT, A
                break

    else:
        fmin1 = fmin2 = fmin

        changedLUT = False

        lut = tag + ".LUT.txt"
        makeLUT(lut, bins, dT, A)

        while True:

            fmin1, new_dT, new_A = optimizeLUT('A', lut, datfile, trials, tag, weight, runtime,
                                                   recursion_level=recursion_level, email=email, verbose=verbose,
                                                   plotenabled=plotenabled)

            if fmin1 < fmin2 and fmin1 < fmin:
                makeLUT(lut, bins, new_dT, new_A)
                changedLUT = True

            fmin2, new_dT, new_A = optimizeLUT('dT', lut, datfile, trials, tag, weight, runtime,
                                               recursion_level=recursion_level, email=email, verbose=verbose,
                                               plotenabled=plotenabled)

            if fmin2 < fmin1 and fmin2 < fmin:
                makeLUT(lut, bins, new_dT, new_A)
                changedLUT = True

            if fmin1 > fmin and fmin2 > fmin:
                if changedLUT:
                    # if it cannot improve it past the fmin, save the best schedule and break the loop
                    break
                else:
                    # if it hasn't improved the given schedule at all, return and break out of the recursion...our job
                    #  is done.
                    return fmin, dT, A
            else:
                fmin = min([fmin1, fmin2])

    if var == "both":
        lut = tag + ".LUT.txt"
        fmin, dT, A = branchLUT(lut, tag, datfile, trials, weight, runtime, recursion_level, email, plotenabled, verbose)
    elif var == 'A':
        A = varmin.copy()
    else:
        dT = varmin.copy()

    if recursion_level == 0:
        if verbose:
            # Print the best updates
            print("Best # updates: " + str(fmin))
        if plotenabled:
            t = np.cumsum(dT)

            t = t - np.ediff1d(t, to_begin=t[0]) / 2.0  # staggers the time so that it falls in between the bins

            plt.plot(t, A)

            plt.ylabel("A-Values")
            plt.xlabel("Time")
            plt.title("A vs. T")
            ax = plt.gca()
            ax.relim()
            ax.autoscale_view()
            plt.draw()
        if email:
            if var == "both":
                msg = "Optimization finished after " + str(datetime.datetime.now() - start) + \
                      ", at " + datetime.datetime.now().strftime("%a %d/%m/%y %H:%M:%S") + \
                      "!\nOptimal dT: {0}\nOptimal A: {1}\nOptimum # updates: {2}\n".format(dTmin, Amin, fmin)
            else:
                msg = "Optimization finished after " + str(datetime.datetime.now() - start) + \
                      ", at " + datetime.datetime.now().strftime("%a %d/%m/%y %H:%M:%S") + \
                      "!\nOptimal " + var + ": " + str(varmin) + "\nOptimum # updates: " + str(fmin) + "\n"
            sendEmail(msg)

    return fmin, dT, A


def getMinimizer(var):
    # Minimize var
    if var == 'dT':
        if xpmt == 1:
            def f(edge, edgeI, tag, filename, trials, dT, A, weight, runtime):
                edges = np.insert(np.cumsum(dT), 0, 0)

                edges[edgeI + 1] = edge

                dT = np.diff(edges)

                return tryLUT(tag, filename, trials, dT, A, weight, runtime, plotenabled, verbose)
        else:
            f = lambda x1, i, x2, a1, a2, a3, a4, a5, a6, v, p: tryLUT(a1, a2, a3, np.insert(x2, i, x1), a4, a5,
                                                                       a6, verbose=v,
                                                                       plotenabled=p)  # rearranging the arguments for dT
    elif var == 'A':
        f = lambda x1, i, x2, a1, a2, a3, a4, a5, a6, v, p: tryLUT(a1, a2, a3, a4, np.insert(x2, i, x1), a5,
                                                                   a6, verbose=v,
                                                                   plotenabled=p)  # rearranging the arguments for A
    elif var == 'both':
        return None
    else:
        raise Exception("Invalid variable argument! Must be \"dT\", \"A\" or \"both\"")
    return f


def branchLUT(lut, tag, datfile, trials, weight, runtime, recursion_level, email, plotenabled, verbose):
    # Get the best var
    bins, dT, A = parseLUT(lut)
    # split dT's largest bin in half
    maxBin = dT.max()
    maxI = dT.argmax()
    dT[maxI] = maxBin / 2.0
    dT = np.insert(dT, maxI, maxBin / 2.0)
    A = np.insert(A, maxI, A[maxI])

    lut = lut.rstrip(".txt") + "." + str(recursion_level) + ".txt"
    makeLUT(lut, bins + 1, dT, A)

    print("Recursing down to level {0}...".format(recursion_level + 1))
    print("dT={0}\nA={1}".format(dT, A))

    return optimizeLUT("both", lut, datfile, trials, tag, weight, runtime, recursion_level=recursion_level + 1,
                              email=email, verbose=verbose, plotenabled=plotenabled)


if __name__ == "__main__":
    sys.exit(main())
