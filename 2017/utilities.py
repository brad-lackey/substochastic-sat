#!/usr/bin/python
import smtplib
import numpy as np


def parseOUT(filename):

    files = []
    times = []
    loops = []
    updates = []

    with open(filename, 'r') as f:
        line = f.readline()
        while len(line) > 0:
            file, _, tStr, loopStr, uStr = line.split()

            files.append(file)
            times.append(float(tStr))
            loops.append(int(loopStr))
            updates.append(int(uStr))

            line = f.readline()

    return files, times, loops, updates


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

"""Returns the percentage of hits and avg runtime as a tuple"""
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