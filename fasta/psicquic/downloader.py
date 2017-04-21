##############################################
# Zhou et al. Cell 2017 -- Parker de Waal 2017
#
# Convert uniprot .txt files to phoscofinder.py compatible sequence files
#
# Usage ex: python uniprot.py uniprot-gpcrs.txt > uniprot-gpcrs.fasta
##############################################

import sys
import re
import os
import numpy as np
import urllib2

inFile = sys.argv[1]
baseUrl = "http://www.uniprot.org/uniprot/"
format = ".txt"

ids = open(inFile, "r")

for id in ids:
	url = baseUrl + str(id).rstrip('\n') + format
	try:
		response = urllib2.urlopen(url)
		html = response.read()
		print html
	except:
		print "failed :", id
