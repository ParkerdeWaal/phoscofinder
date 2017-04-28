##############################################
# Zhou et al. Cell 2017 -- Parker de Waal 2017
#
# Convert uniprot .txt files to phoscofinder.py compatible sequence files
#
# Usage ex: python uniprot.py Type(GPCR|non) location(cTail|ICL3) uniprot-gpcrs.txt > uniprot-gpcrs.fasta
##############################################

import sys
import re
import os
import numpy as np

type = str(sys.argv[1])
loc = str(sys.argv[2])
inFile = sys.argv[3]

sequences = open(inFile, "r")

idFound = 0
seqIn = 0
isGPCR = "non"


cytoFormat = np.dtype([('FT','S8'),('FTT','S12'),('begin', 'i'),('end','i')])
cytoRanges = np.empty([0,3],dtype=cytoFormat)


def processCyto(ID,genesequence,cytoRanges):
	global loc
	if cytoRanges[-1][0] == "TOPO_DOM" and cytoRanges[-1][1] == "Cytoplasmic." and "TRANSMEM" in cytoRanges['FT']:
		if loc == "cTail":
			begin = cytoRanges[-1][2]
			end = cytoRanges[-1][3]
		elif loc == "ICL3":
			begin = cytoRanges[-5][2]
			end = cytoRanges[-5][3]
		print ID,begin,end,sequence[begin-1:end]

	
for line in sequences:
	# determine GPCR status
        if "GPCR_" in line:
                isGPCR = "GPCR"

	if line[0:2] == "KW" and "G-protein coupled receptor" in line:
		isGPCR = "GPCR"

	# Protein is over, call process
        if line[0:2] == "//":
                seqIn = 0
                idFound = 0
		if cytoRanges.size > 0 and isGPCR == type:
			processCyto(ID,sequence,cytoRanges)
		isGPCR = "non"
		cytoRanges = np.empty([0,3],dtype=cytoFormat)
	# protein beings
	if line[0:2] == "ID" and idFound == 0:
		idFound = 1
		ID = line.split()[1]
	# extract domains		
	if line[0:13] == "FT   TRANSMEM" or line[0:13] == "FT   TOPO_DOM":
		if str(line.split()[2]) == "?":
			continue	
		FT = str(line.split()[1])
		FTT = str(line.split()[4])
		tempBegin = int(line[15:20])
		tempEnd   = int(line[21:30])
		cytoRanges = np.append(cytoRanges, np.array([(FT,FTT,tempBegin,tempEnd)], dtype=cytoFormat))
	
	# aggregate sequence data
	if seqIn == 1:	
		seqLine = str(line.replace(' ','').replace('\n',''))
		sequence+=seqLine
	# check for sequence begin
        if line[0:2] == "SQ":
                seqIn = 1
		sequence = ""




