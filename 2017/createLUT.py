#! /usr/bin/python
# Make a tuple out of a string
def make_tuple(s):
	return tuple(s[1:-1].split(','))


print("Enter Filename:")
filename = raw_input()
print("Enter time spacing:")
dT = int(raw_input())
print("Enter range of values as (val0, valN):")
val0, valN = make_tuple(raw_input())
val0, valN = float(val0), float(valN)
print("Enter # Bins:")
bins = int(raw_input())

deltaV = (valN-val0)/float(bins)

with open(filename, 'w') as f:
	f.write(str(bins+1) + '\n')
	for i in range(bins+1):
		f.write(str(dT) + '\t' + str(val0 + deltaV*i) + '\n')

print("LUT complete. Written to {0}".format(filename))