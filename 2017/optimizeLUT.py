#!/usr/bin/python
import sys
from createLUT import makeLUT
from subprocess import check_call
import numpy as np
from scipy.optimize import basinhopping
import matplotlib.pyplot as plt
import smtplib

num_iter = 0
email = False
verbose = False

"""Returns the percentage of hits and avg runtime as a tuple"""
def parseTXT(txtfile):
    last = ''
    with open(txtfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            last = line
            line = f.readline()

    # Parse last line
    _, hitStr, _, tStr = last.split()

    fraction = map(float, hitStr[0:hitStr.rfind('(')].split('/'))
    hit = fraction[0]/fraction[1]

    t = float(tStr.rstrip("s"))

    return hit, t

"""Returns the avg time of a set of conf files using given LUT"""
def tryLUT(tag, filename, trials, dT, A, weight, runtime, verbose):
    if len(dT) != len(A):
        raise Exception("Vectors dT and A are not the same length!")

    bins = len(dT)
    lut = tag + ".LUT." + str(bins) + ".txt"

    makeLUT(lut, bins, dT, A)

    args = []
    args.append('./testrun.pl')  # the program to run
    args.append('./ssmc')
    args.append(lut)
    args.append(filename)
    args.append(trials)
    args.append(tag)

    if weight and runtime:
        args.append(weight)
        args.append(runtime)

    # returns 0 if successful otherwise throws error
    check_call(args)

    txtfile = tag + ".txt"
    hits, avg_time = parseTXT(txtfile)

    penalty = 10000
    if hits < 1:
        avg_time = penalty

    if verbose:
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

        global num_iter
        num_iter += 1

        print("# Iterations: " + str(num_iter))
        print("Tried dT=" + str(dT) + ", A=" + str(A) + " with a time of " + str(avg_time))

    return avg_time


def sendEmail(msg):
    server = smtplib.SMTP("smtp.gmail.com", 587)
    server.starttls()
    server.login("email.notifier.bryanluu@gmail.com", "7788382652")
    server.sendmail("email.notifier.bryanluu@gmail.com", "bryanluu30794@gmail.com", msg)
    server.quit()


def main():
    args = sys.argv
    if '-m' in args:
        global email
        email = True
        args.remove('-m')
    if '-v' in args:
        global verbose
        verbose = True
        args.remove('-v')
    if len(args) == 5 or len(args) == 7:
        filename = args[1]
        bins = int(args[2])
        trials = args[3]
        tag = args[4]
        weight = None
        runtime = None

        if len(args) == 7:
            weight = args[5]
            runtime = args[6]

    else:
        print("Usage: ./optimizeLUT [-v] [-m] <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
        return 1

    if verbose:
        # Turn on interactive plotting
        plt.ion()

    # Guesses for dT and A
    dT = np.ones(bins)
    A = np.linspace(1, 0, bins)

    # Minimize A
    bnds = tuple((0,1) for x in A)  # A values lie between (0,1)
    f_A = lambda x, a1, a2, a3, a4, a5, a6, a7: tryLUT(a1, a2, a3, a4, x, a5, a6, a7)  # rearranging the arguments for A

    def printf(x, f, accept):
        if verbose:
            print("#################### Found minimum: " + str(f) + " ####################")
            global num_iter
            num_iter = 0
        if email:
            msg = "Found minimum: " + str(f)
            sendEmail(msg)

    res = basinhopping(f_A, A, stepsize=0.1, callback=printf,
                    minimizer_kwargs={'method':'L-BFGS-B', 'options':{'eps':0.005, 'ftol':0.001, 'maxfun':bins*10}, 'bounds':bnds,
                                      'args':(tag, filename, trials, dT, weight, runtime, verbose)})

    # Store the best A
    A = res.x
    lut = tag + ".OPTIMAL.A.txt"
    makeLUT(lut, bins, dT, A)

    if verbose:
        # Print the best time
        print("Best time: " + str(res.fun))

        ax = plt.gca()
        t = np.cumsum(dT)
        t = t - np.ediff1d(t, to_begin=t[0])/2.0  # staggers the time so that it falls in between the bins
        ax.plot(t, A)
        ax.ylabel("A-Values")
        ax.xlabel("Time")
        ax.title("A vs. T")
        ax.relim()
        ax.autoscale_view()
        plt.draw()

    if email:
        msg = "Optimization finished!\nOptimal A: " + str(A) + "\nOptimum value: " + str(res.fun) + "\n"
        sendEmail(msg)


    return 0

if __name__ == "__main__":
    sys.exit(main())