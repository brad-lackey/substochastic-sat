#!/usr/bin/python

from optimizeLUT import parseLUT, sendEmail, getABounds, plotLUT, plotPsize, tryLUT
from scipy import stats
from subprocess32 import check_call, TimeoutExpired
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fminbound
import datetime
import sys
from createLUT import makeLUT
from histAnalysis import parseOUT
from simanneal import Annealer

BOUND_CAP = 0.1  # cap on the bounds
BOUND_MULTIPLIER = 1.1  # fraction over which the bound can extend
UPDATE_PENALTY = 10000000  # penalty to give scripts which timeout
N_ITERS_CAP = 5  # max number of optimization iterations
RECURSION_LIMIT = 5  # max levels optimizer can branch LUT
THRESHOLD = 0.25  # min threshold before accepting new minimum


class Optimizer(Annealer):
    """Optimize using simulated annealing"""

    def __init__(self, var, row, state, other1, other2, tag, datfile, trials, plotenabled=False, verbose=False):
        self.other1, self.other2 = other1, other2
        self.tag = tag
        self.datfile = datfile
        self.trials = trials
        self.var = var
        self.row = row
        self.forward = True
        self.plotenabled, self.verbose = plotenabled, verbose
        super(Optimizer, self).__init__(state)

    def move(self):
        """Perturb the current index randomly"""
        bins = len(self.state)

        if self.var == "A":
            lbound, ubound = 0.1, 1.0
        elif self.var == "dT":
            lbound, ubound = 0.1, 100.0
        else:
            # if var == psize
            lbound, ubound = max([self.state[self.row]-4, 8]), self.state[self.row]+4

        scale = 0.1
        diff = (ubound-lbound)*0.5

        # perturb the given row
        self.state[self.row] += np.random.randn() * scale * diff

        # check bounds
        self.state[self.row] = min([ubound, self.state[self.row]])
        self.state[self.row] = max([lbound, self.state[self.row]])



    def energy(self):
        e = 10000000000

        if self.var == "A":
            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.other1, self.state, self.other2,
                       plotenabled=self.plotenabled, verbose=self.verbose)
        elif self.var == "dT":
            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.state, self.other1, self.other2,
                       plotenabled=self.plotenabled, verbose=self.verbose)
        elif self.var == "psize":
            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.other1, self.other2, self.state,
                       plotenabled=self.plotenabled, verbose=self.verbose)

        return e


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
    fmin = tryLUT(var, tag, datfile, trials, dT, A, psize, weight, runtime, plotenabled, verbose)

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
            fval = 0.0

            if xpmt != 0:
                np.random.shuffle(indices)  # shuffles the indices array for random choice of index to optimize

            minFound = False


            for row in indices:
                if var == "dT":
                    varvector, other1, other2 = dT, A, psize
                elif var == "A":
                    varvector, other1, other2 = A, dT, psize
                else:
                    varvector, other1, other2 = psize, dT, A

                opt = Optimizer(var, row, varvector.tolist(), other1, other2, tag, datfile, trials, plotenabled, verbose)

                opt.copy_strategy = "slice"

                opt.Tmax = 5000000/float(trials)  # Max (starting) temperature
                opt.Tmin = 1000000/float(trials)      # Min (ending) temperature
                opt.steps = 100   # Number of iterations
                opt.updates = 100   # Number of updates (by default an update prints to stdout)

                vlist, fval = opt.anneal()

                varvector = np.array(vlist)

                x0 = varvector[row]

                if fval < fmin:
                    fmin = fval

                    _, _, _, best_updates = parseOUT(tag + ".out")

                    varmin = varvector.copy()

                    printf(fmin)

                    # Store the best var
                    lut = tag + ".OPTIMAL." + var + ".lut"
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
                        "---------- Found {0}[{1}]={2}".format(var, row, x0) + " at updates " + str(fval) +
                        ", {0}/{1} iterations ----------".format(i + 1, N_ITERS_CAP))

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

        lut = tag + ".lut"
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
        lut = tag + ".lut"
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

    lut = lut.rstrip(".lut") + "." + str(recursion_level) + ".lut"

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
