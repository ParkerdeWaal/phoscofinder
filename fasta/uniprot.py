##############################################
# Zhou et al. Cell 2017 -- Parker de Waal 2017
#
# Convert uniprot .txt files to phoscofinder.py compatible .fasta files
#
# Usage ex: python uniprot.py uniprot-gpcrs.txt > uniprot-gpcrs.fasta
##############################################

import sys
import re
import os
import numpy as np
import urllib2

inFile = sys.argv[1]

sequences = open(inFile, "r")

idFound = 0
seqIn = 0


cytoFormat = np.dtype([('FT','S8'),('FTT','S12'),('begin', 'i'),('end','i')])
cytoRanges = np.empty([0,3],dtype=cytoFormat)


def processCyto(ID,genesequence,cytoRanges):
	if cytoRanges[-1][0] == "TOPO_DOM" and cytoRanges[-1][1] == "Cytoplasmic.":
		begin = cytoRanges[-1][2]
		end = cytoRanges[-1][3]
		print ID,begin,end,sequence[begin-1:end]

	
for line in sequences:
	## Protein is over
        if line[0:2] == "//":
                seqIn = 0
                idFound = 0
		if cytoRanges.size > 0:
			processCyto(ID,sequence,cytoRanges)
		cytoRanges = np.empty([0,3],dtype=cytoFormat)

	if line[0:2] == "ID" and idFound == 0:
		idFound = 1
		ID = line.split()[1]		
	if line[0:13] == "FT   TRANSMEM" or line[0:13] == "FT   TOPO_DOM":
		if str(line.split()[2]) == "?":
			continue
		
		FT = str(line.split()[1])
		FTT = str(line.split()[4])
		tempBegin = int(line[15:20])
		tempEnd   = int(line[21:30])
		cytoRanges = np.append(cytoRanges, np.array([(FT,FTT,tempBegin,tempEnd)], dtype=cytoFormat))

	if seqIn == 1:	
		seqLine = str(line.replace(' ','').replace('\n',''))
		sequence+=seqLine

        if line[0:2] == "SQ":
                seqIn = 1
		sequence = ""




