#!/usr/bin/python
import sys
from createLUT import makeLUT
from subprocess import check_call
import numpy as np
from scipy.optimize import fminbound
import matplotlib.pyplot as plt
import smtplib

num_iter = 0
email = False
verbose = False

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

    return dT, A

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
def tryLUT(tag, filename, trials, dT, A, weight, runtime):
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
    if len(args) == 7 or len(args) == 9:
        lutfile = args[1]
        var = args[2]
        datfile = args[3]
        bins = int(args[4])
        trials = args[5]
        tag = args[6]
        weight = None
        runtime = None

        if len(args) == 9:
            weight = args[7]
            runtime = args[8]

    else:
        print("Usage: ./optimizeLUT dT|A [-v] [-m] <initialLUT> <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n")
        return 1

    if verbose:
        # Turn on interactive plotting
        plt.ion()

    # Load initial conditions for dT and A
    dT, A = parseLUT(lutfile)

    # Minimize var
    if var == 'dT':
        f_A = lambda x1, i, x2, a1, a2, a3, a4, a5, a6: tryLUT(a1, a2, a3, np.insert(x2, i, x1), a4, a5, a6)  # rearranging the arguments for A
    elif var == 'A':
        f_A = lambda x1, i, x2, a1, a2, a3, a4, a5, a6: tryLUT(a1, a2, a3, a4, np.insert(x2, i, x1), a5, a6)  # rearranging the arguments for A
    else:
        raise Exception("Invalid variable argument! Must be \"dT\" or \"A\"")

    def printf(f):
        if verbose:
            print("#################### Found new minimum: " + str(f) + " ####################")
            global num_iter
            num_iter = 0
        if email:
            msg = "Found new minimum: " + str(f)
            sendEmail(msg)

    # res = basinhopping(f_A, A, stepsize=0.1, callback=printf,
    #                 minimizer_kwargs={'method':'TNC', 'options':{'eps':0.005, 'ftol':0.01}, 'bounds':bnds,
    #                                   'args':(tag, filename, trials, dT, weight, runtime)})

    # Loop through each index of dT and minimize accordingly
    N = 10
    fmin = 1000
    varmin = -1
    for n in range(N):
        fval = 0
        for i in range(bins):
            if var == "dT":
                x0, fval, ierr, numfunc = fminbound(f_A, 0, 1, args=(i, np.delete(dT,i), tag, datfile, trials, A, weight, runtime),
                      full_output=True, xtol=0.01)
                dT[i] = x0
            else:
                x0, fval, ierr, numfunc = fminbound(f_A, 0, 1, args=(i, np.delete(A,i), tag, datfile, trials, dT, weight, runtime),
                                                    full_output=True, xtol=0.01)
                A[i] = x0

            if fval < fmin:
                fmin = fval

                if var == "dT":
                    varmin = dT.copy()
                else:
                    varmin = A.copy()

                printf(fmin)

                # Store the best var
                lut = tag + ".OPTIMAL." + var + ".txt"
                if var == "dT":
                    makeLUT(lut, bins, varmin, A)
                else:
                    makeLUT(lut, bins, dT, varmin)

                if verbose:
                    plt.savefig(tag + ".OPTIMAL." + var + ".png")

            if verbose:
                print("---------- Found " + str(x0) + " at " + str(fval) + " after " + str(numfunc) + " iterations ----------")


    if verbose:
        # Print the best time
        print("Best time: " + str(fmin))

        ax = plt.gca()

        if var == 'dT':
            t = np.cumsum(varmin)
        else:
            t = np.cumsum(dT)

        t = t - np.ediff1d(t, to_begin=t[0])/2.0  # staggers the time so that it falls in between the bins

        if var == 'A':
            plt.plot(t, varmin)
        else:
            plt.plot(t, A)

        plt.ylabel("A-Values")
        plt.xlabel("Time")
        plt.title("A vs. T")
        plt.relim()
        plt.autoscale_view()
        plt.draw()

    if email:
        msg = "Optimization finished!\nOptimal " + var + ": " + str(varmin) + "\nOptimum value: " + str(fmin) + "\n"
        sendEmail(msg)


    return 0

if __name__ == "__main__":
    sys.exit(main())