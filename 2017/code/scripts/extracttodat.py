from bs4 import BeautifulSoup
import urllib
import re
import unicodedata
import os
import fnmatch

r = urllib.urlopen('http://maxsat.ia.udl.cat/detailed/incomplete-ms-crafted-table.html').read()
soup = BeautifulSoup(r,"lxml")
rows = soup.findAll('tr')
for row in rows:
	files = row.findChildren('td', {"bgcolor" : "#DFDFFF"})
	solutions = row.findChildren('td', {"bgcolor" : "#FFFFFF"})
	f_name = ""
	for file in files:
		value = file.string
#		print value
		for root, dirnames, filenames in os.walk('./ms_crafted'):
                	for filename in fnmatch.filter(filenames, value):
                        	f_name = os.path.join(root,filename)
	for solution in solutions:
#		print solution
                value = solution.findChildren("font")
#		print value
		sol = solution.font.contents[0].string.encode('ascii','ignore')
		time = "T = 3.00"
#		if len(solution.font.contents) > 1:
#			time = solution.font.contents[1].string.encode('ascii','ignore')
		print f_name + "   " + sol + "   " + time
#	break
		
#	for root, dirnames, filenames in os.walk('.'):
#		for filename in fnmatch.filter(filenames, '*.wcnf'):
#			print os.path.join(root,filename)

#        for root, dirnames, filenames in os.walk('.'):
#                for filename in fnmatch.filter(filenames, '*.cnf'):
#                        print os.path.join(root,filename)

