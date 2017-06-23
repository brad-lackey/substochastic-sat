#!/usr/bin/python
import smtplib
import numpy as np

"""Returns the ratio of a given CNF file"""
def parseCNF(cnf):
    with open(cnf, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            if line.startswith('p'):

                if cnf.endswith(".cnf"):
                    p, fmt, varStr, cStr = line.split()
                elif cnf.endswith(".wcnf"):
                    p, fmt, varStr, cStr, pStr = line.split()
                else:
                    raise Exception("Invalid CNF file: {0}".format(cnf))

                var = float(varStr)
                clauses = float(cStr)

                ratio = (var/clauses)

                return ratio

            line = f.readline()

        return None



"""Returns the files, optima and times, in the given DAT file"""
def parseDAT(datfile):

    files = []
    optima = []
    times = []

    with open(datfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            filename, oLbl, eq, oStr, tLbl, eq, tStr = line.split()

            files.append(filename)
            optima.append(int(oStr))
            times.append(float(tStr))

            line = f.readline()

    return files, optima, times


"""Parse a .out file from testrun, returning the files, times, loops and updates."""
def parseOUT(filename):

    files = []
    optima = []
    times = []
    loops = []
    updates = []

    with open(filename, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            filename, oStr, tStr, loopStr, uStr = line.split()

            files.append(filename)
            optima.append(int(oStr))
            times.append(float(tStr))
            loops.append(int(loopStr))
            updates.append(int(uStr))

            line = f.readline()

    return files, optima, times, loops, updates


"""Returns the dT and A vectors from a LUT file as a tuple"""
def parseLUT(lutfile):

    with open(lutfile, 'r') as f:
        line = f.readline()
        bins = int(line)
        dT = np.zeros(bins)
        A = np.zeros(bins)
        psize = np.zeros(bins)
        row = 0
        while len(line) > 0:
            line = f.readline()
            if len(line) > 0:
                dT[row], A[row], psize[row] = line.rstrip('\n').split('\t')
                row += 1
        if row != bins:
            raise Exception("Invalid LUT file format!")

    return bins, dT, A, psize


"""Returns the percentage of hits, avg runtime, and factor as a tuple"""
def parseTXT(txtfile):
    last = ''
    with open(txtfile, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            last = line
            line = f.readline()

    # Parse last line
    _, hitStr, _, tStr, lStr, _, uStr, _, fStr, _ = last.split()

    fraction = map(float, hitStr[0:hitStr.rfind('(')].split('/'))
    hit = fraction[0]/fraction[1]
    loops = float(lStr)
    t = float(tStr.rstrip("s"))
    updates = float(uStr)
    factor = float(fStr)

    return hit, updates, factor


def sendEmail(msg):
    server = smtplib.SMTP("smtp.gmail.com", 587)
    server.starttls()
    server.login("email.notifier.bryanluu@gmail.com", "7788382652")
    server.sendmail("email.notifier.bryanluu@gmail.com", "bryanluu30794@gmail.com", '\n'+msg)
    server.quit()