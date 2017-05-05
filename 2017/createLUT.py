#! /usr/bin/python
# Make a tuple out of a string
def make_tuple(s):
	return tuple(s[1:-1].split(','))


print("Enter Filename:")
filename = raw_input()
print("Enter Beginning Point as (time, val):")
time0, val0 = make_tuple(raw_input())
time0, val0 = int(time0), float(val0)
print("Enter End Point as (time, val):")
timeN, valN = make_tuple(raw_input())
timeN, valN = int(timeN), float(valN)
print("Enter # Bins:")
bins = int(raw_input())

deltaT = (timeN-time0)/float(bins)
deltaV = (valN-val0)/float(bins)

with open(filename, 'w') as f:
	f.write(str(bins+1) + '\n')
	for i in range(bins+1):
		f.write(str(int(time0 + deltaT*i)) + '\t' + str(val0 + deltaV*i) + '\n')

print("LUT complete. Written to {0}".format(filename))