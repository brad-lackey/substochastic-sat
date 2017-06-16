#!/usr/bin/python

from optimizeLUT import parseLUT, parseTXT, sendEmail, plotLUT, plotPsize
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
import math

BOUND_CAP = 0.1  # cap on the bounds
BOUND_MULTIPLIER = 1.1  # fraction over which the bound can extend
UPDATE_PENALTY = 10000000  # penalty to give scripts which timeout
N_ITERS_CAP = 1  # max number of optimization iterations
RECURSION_LIMIT = 5  # max levels optimizer can branch LUT
THRESHOLD = 0.25  # min threshold before accepting new minimum


class Optimizer(Annealer):
    """Optimize using simulated annealing"""

    def __init__(self, var, state, other1, other2, tag, datfile, trials, plotenabled=False, verbose=False):
        self.other1, self.other2 = other1, other2
        self.tag = tag
        self.datfile = datfile
        self.trials = trials
        self.var = var
        self.forward = True
        self.plotenabled, self.verbose = plotenabled, verbose

        if var == "all":
            fullstate = state + other1 + other2
        else:
            fullstate = state
        super(Optimizer, self).__init__(fullstate)

    def move(self):
        """Perturb the current index randomly"""
        bins = len(self.other1)

        for row in range(bins):
            if self.var == "A":
                lbound, ubound = 0.1, 1.0
                mu_diff, sigma = 0.05, 0.05
            elif self.var == "dT":
                lbound, ubound = 0.1, 100.0
                mu_diff, sigma = 0.5, 0.5
            elif self.var == "psize":
                lbound, ubound = 16, 128
                mu_diff, sigma = 5, 5
            else:
                self.walk(row, 0.5, 0.5, 0.1, 100.0)  # perturb dT
                self.walk(bins+row, 0.05, 0.05, 0.1, 1.0)  # perturb A
                self.walk(bins+bins+row, 5, 5, 16, 128)
                continue

            self.walk(row, mu_diff, sigma, lbound, ubound)

    def walk(self, row, mu_diff, sigma, lbound, ubound):
        # perturb the given row
        coin = np.random.randint(2)
        mu = self.state[row]
        if coin == 0:
            self.state[row] = sigma * np.random.randn() + mu - mu_diff
        else:
            self.state[row] = sigma * np.random.randn() + mu + mu_diff

        # check bounds
        self.state[row] = min([ubound, self.state[row]])
        self.state[row] = max([lbound, self.state[row]])

    def energy(self):
        if self.var == "A":
            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.other1, self.state, self.other2,
                       plotenabled=self.plotenabled, verbose=self.verbose)
        elif self.var == "dT":
            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.state, self.other1, self.other2,
                       plotenabled=self.plotenabled, verbose=self.verbose)
        elif self.var == "psize":
            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.other1, self.other2, self.state,
                       plotenabled=self.plotenabled, verbose=self.verbose)
        else:
            bins = len(self.other1)

            e = tryLUT(self.var, self.tag, self.datfile, self.trials, self.state[:bins],
                       self.state[bins:(bins+bins)], self.state[(bins+bins):],
                       plotenabled=self.plotenabled, verbose=self.verbose)

        return e


"""Returns the factor of a set of conf files using given LUT"""
def tryLUT(var, tag, filename, trials, dT, A, psize, weight=None, runtime=None, plotenabled=False, verbose=False):
    if len(dT) != len(A) or len(psize) != len(dT) or len(psize) != len(A):
        raise Exception("Vectors dT, A and psize are not the same length!")

    bins = len(dT)
    lut = tag + ".lut"

    psize = map(round, psize)
    psize = map(int, psize)

    makeLUT(lut, bins, dT, A, psize)

    args = []
    args.append("./testrun.pl")  # the program to run
    args.append("./ssmc")
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
    hits, updates, factor = parseTXT(txtfile)

    if hits < 1:
        updates += factor*UPDATE_PENALTY

    if plotenabled:
        if var == 'psize':
            plotPsize(dT, psize)
        else:
            plotLUT(dT, A)

    if verbose:
        print("Tried dT=" + str(dT) + ", A=" + str(A) + ", Psize=" + str(psize) + "  with factor=" + str(factor))

    return factor


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
    if len(args) == 6 or len(args) == 8:
        global var
        var = args[1]

        lutfile = args[2]
        datfile = args[3]
        trials = args[4]
        tag = args[5]
        weight = None
        runtime = None

        if len(args) == 8:
            weight = args[6]
            runtime = args[7]

    else:
        print("Usage: ./annealer.py dT|A|psize|all [-v] [-m] [-p] <initialLUT> <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
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

    def printf(fval):
        if verbose:
            print("#################### Found new minimum: " + str(fval) + " ####################")
        if email:
            msg = "Found new minimum: " + str(fval)
            sendEmail(msg)

    # set initial minimum to initial LUT performance
    fmin = tryLUT(var, tag, datfile, trials, dT, A, psize, weight, runtime, plotenabled, verbose)

    if var == 'dT':
        varmin = dT.copy()
    elif var == 'A':
        varmin = A.copy()
    elif var == "psize":
        varmin = psize.copy()
    elif var == "all":
        pass
    else:
        raise Exception("Invalid variable argument! Must be \"dT\", \"A\", \"psize\" or \"all\"")

    if var != "all":

        if var == "dT":
            varvector, other1, other2 = dT, A, psize
        elif var == "A":
            varvector, other1, other2 = A, dT, psize
        else:
            varvector, other1, other2 = psize, dT, A

        opt = Optimizer(var, varvector.tolist(), other1, other2, tag, datfile, trials, plotenabled, verbose)

        # the state vector can simply be copied by slicing
        opt.copy_strategy = "slice"

        opt.Tmax = 10  # Max (starting) temperature
        opt.Tmin = 0.1      # Min (ending) temperature
        opt.steps = 1000   # Number of iterations

        try:
            vlist, fval = opt.anneal()
        except Exception:
            vlist, fval = opt.best_state, opt.best_energy

        varvector[:] = vlist[:]



    else:
        varvector, other1, other2 = dT, A, psize

        opt = Optimizer(var, varvector.tolist(), other1.tolist(), other2.tolist(), tag, datfile, trials, plotenabled, verbose)

        # the state vector can simply be copied by slicing
        opt.copy_strategy = "slice"

        opt.Tmax = 10  # Max (starting) temperature
        opt.Tmin = 0.1      # Min (ending) temperature
        opt.steps = 1000   # Number of iterations

        try:
            vlist, fval = opt.anneal()
        except Exception:
            vlist, fval = opt.best_state, opt.best_energy

        varvector[:] = vlist[:bins]
        other1[:] = vlist[bins:(bins+bins)]
        other2[:] = vlist[(bins+bins):]


    if fval < fmin:
        fmin = fval

        varmin = varvector.copy()

        printf(fmin)

        # Store the best var
        lut = tag + ".OPTIMAL." + var + ".lut"
        if var == "dT":
            makeLUT(lut, bins, varmin, A, psize)
        elif var == "psize":
            makeLUT(lut, bins, dT, A, varmin)
        elif var == "A":
            makeLUT(lut, bins, dT, varmin, psize)
        else:
            makeLUT(lut, bins, dT, A, psize)

        if plotenabled:
            plt.savefig(tag + ".OPTIMAL." + var + ".png")


    if var == 'A':
        A = varmin.copy()
    elif var == 'psize':
        psize = varmin.copy()
    else:
        dT = varmin.copy()

    if verbose:
        # Print the best factor
        print("Best factor: " + str(fmin))
    if plotenabled:
        if var == "psize":
            plotPsize(dT, psize)
        else:
            plotLUT(dT, A)
    if email:
        if var == "both":
            msg = "Optimization finished after " + str(datetime.datetime.now() - start) + \
                  ", at " + datetime.datetime.now().strftime("%a %d/%m/%y %H:%M:%S") + \
                  "!\nOptimal dT: {0}\nOptimal A: {1}\nOptimum factor: {2}\n".format(dT, A, fmin)
        else:
            msg = "Optimization finished after " + str(datetime.datetime.now() - start) + \
                  ", at " + datetime.datetime.now().strftime("%a %d/%m/%y %H:%M:%S") + \
                  "!\nOptimal " + var + ": " + str(varmin) + "\nOptimum factor: " + str(fmin) + "\n"
        sendEmail(msg)

    return fmin, dT, A, psize


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
