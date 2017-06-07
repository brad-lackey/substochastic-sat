#!/usr/bin/python
import sys

# Make a tuple out of a string
def make_tuple(s, d_type):
    return map(d_type, tuple(s[1:-1].split(',')))

# create the LUT table
def makeLUT(filename, bins, dT, A, psize):
    dT = map(str, dT)
    A = map(str, A)
    psize = map(str, psize)

    with open(filename, 'w') as f:
        f.write(str(bins) + '\n')
        for i in range(bins):
            f.write(dT[i] + '\t' + A[i] + '\t' + psize[i] + '\n')

def main():

    # Whether the script accepts user input or command line input
    cmd_line = len(sys.argv) > 1

    if cmd_line:
        filename = sys.argv[1]
    else:
        print("Enter Filename:")
        filename = raw_input().rstrip()

    if cmd_line:
        mode = int(sys.argv[2])
    else:
        print("Enter input mode (0:uniform interpolation, 1:manual):")
        mode = int(raw_input())

    # Mode 0 means uniform spacing, 1 is manual entry of dT and A-values
    if mode is 0:
        if cmd_line:
            dT = sys.argv[3]
        else:
            print("Enter time spacing:")
            dT = raw_input()

        # float verification
        float(dT)

        if cmd_line:
            val0, valN = map(float, sys.argv[4:6])
        else:
            print("Enter range of values as (val0, valN):")
            val0, valN = make_tuple(raw_input(), float) # keep val as a float for calculation

        if cmd_line:
            psize = int(sys.argv[6])
        else:
            print("Enter the uniform population size:")
            psize = int(raw_input())

        if cmd_line:
            bins = int(sys.argv[7])
        else:
            print("Enter # Bins:")
            bins = int(raw_input())

        deltaA = (valN-val0)/(bins-1)

        A = [val0]
        dT = [dT]
        psize = [psize]
        for i in range(bins - 1):
            dT.append(dT[0])
            A.append(val0 + (i+1)*deltaA)
            psize.append(psize[0])

    elif mode is 1:
        if cmd_line:
            bins = int(sys.argv[3])
        else:
            print("Enter # Bins:")
            bins = int(raw_input())

        input_valid = False

        dT = []
        while not input_valid:
            if cmd_line:
                dT = sys.argv[4:4+bins]
            else:
                print("Enter time spacing (separated by spaces, if single, then uniform):")
                dT = raw_input().split()

            # float verification
            map(float, dT)

            if len(dT) == 1:
                # copy first value to the whole list
                for i in range(bins):
                    dT.append(dT[0])
                input_valid = True
            elif len(dT) == bins:
                # n-th value is n-th dT
                input_valid = True
            else:
                if cmd_line:
                    raise Exception("Invalid number of entries!")
                else:
                    input_valid = False
                    print("Invalid number of entries!")

        input_valid = False

        A = []
        while not input_valid:
            if cmd_line:
                A = sys.argv[4+bins: 4+bins+bins]
            else:
                print("Enter A-values (separated by spaces, if single, then uniform):")
                A = raw_input().split()

            # float verification
            map(float, A)

            if len(A) == 1:
                # copy first value to the whole list
                for i in range(bins):
                    A.append(A[0])
                input_valid = True
            elif len(A) == bins:
                # n-th value is n-th dT
                input_valid = True
            else:
                if cmd_line:
                    raise Exception("Invalid number of entries!")
                else:
                    input_valid = False
                    print("Invalid number of entries!")

        input_valid = False

        psize = []
        while not input_valid:
            if cmd_line:
                psize = sys.argv[4+bins+bins:4+bins+bins+bins]
            else:
                print("Enter population sizes (separated by spaces, if single, then uniform):")
                psize = raw_input().split()

            # float verification
            map(int, psize)

            if len(psize) == 1:
                # copy first value to the whole list
                for i in range(bins):
                    psize.append(psize[0])
                input_valid = True
            elif len(psize) == bins:
                # n-th value is n-th dT
                input_valid = True
            else:
                if cmd_line:
                    raise Exception("Invalid number of entries!")
                else:
                    input_valid = False
                    print("Invalid number of entries!")

    else:
        raise Exception("Mode is not 0 or 1, invalid!")

    makeLUT(filename, bins, dT, A, psize)

    print("LUT complete. Written to {0}".format(filename))

    return 0

if __name__ == "__main__":
    sys.exit(main())