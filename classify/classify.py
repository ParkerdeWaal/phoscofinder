##############################################
# Zhou et al. Cell 2017 -- Parker de Waal 2017
#
# Convert uniprot .txt files to phoscofinder.py compatible sequence files
#
# Usage ex: python uniprot.py Type(GPCR|non) uniprot-gpcrs.txt > uniprot-gpcrs.fasta
##############################################

import sys
import re
import csv
import pickle
import os
import numpy as np

type = str(sys.argv[1])
#matchFile = sys.argv[2]
inFile = sys.argv[2]

sequences = open(inFile, "r")

idFound = 0
seqIn = 0
isGPCR = "non"


cytoFormat = np.dtype([('FT','S8'),('FTT','S12'),('begin', 'i'),('end','i')])
rangesFormat = np.dtype([('gene','S11'),('icl3', 'f'),('ctail','f')])
cytoRanges = np.empty([0,3],dtype=cytoFormat)
sizeRanges = np.empty([0,3],dtype=rangesFormat)


def processCyto(ID,genesequence,cytoRanges):
	global sizeRanges
	if cytoRanges[-1][0] == "TOPO_DOM" and cytoRanges[-1][1] == "Cytoplasmic." and "TRANSMEM" in cytoRanges['FT']:
		ICL3l = cytoRanges[-5][3] - cytoRanges[-5][2] +1
		CtailL = cytoRanges[-1][3] - cytoRanges[-1][2] +1
		print ID,ICL3l,CtailL
		sizeRanges = np.append(sizeRanges, np.array([(str(ID),ICL3l,CtailL)], dtype=rangesFormat))
	
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
		if cytoRanges.size > 0 and isGPCR == type: #and ID in matchList:
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



for i in range(0,sizeRanges.size):
    sizeRanges['icl3'][i] = np.log(sizeRanges['icl3'][i])
    sizeRanges['ctail'][i] = np.log(sizeRanges['ctail'][i])

longLoop = []
for i in range(0,sizeRanges.size):
    if int(sizeRanges['icl3'][i]) > np.mean(sizeRanges['icl3']) + np.std(sizeRanges['icl3']) and int(sizeRanges['ctail'][i]) < np.mean(sizeRanges['ctail']) + np.std(sizeRanges['ctail']):
        longLoop.append(sizeRanges['gene'][i])


longTail = []
for i in range(0,sizeRanges.size):
    if int(sizeRanges['ctail'][i]) > np.mean(sizeRanges['ctail']) + np.std(sizeRanges['ctail']) and int(sizeRanges['icl3'][i]) < np.mean(sizeRanges['icl3']) + np.std(sizeRanges['icl3']):
        longTail.append(sizeRanges['gene'][i])

longLoopTail = []
for i in range(0,sizeRanges.size):
    if int(sizeRanges['ctail'][i]) > np.mean(sizeRanges['ctail']) + np.std(sizeRanges['ctail']) and int(sizeRanges['icl3'][i]) > np.mean(sizeRanges['icl3']) + np.std(sizeRanges['icl3']):
        longLoopTail.append(sizeRanges['gene'][i])

normalLoopTail = []
for i in range(0,sizeRanges.size):
    if sizeRanges['gene'][i] not in longLoop and sizeRanges['gene'][i] not in longTail and sizeRanges['gene'][i] not in longLoopTail:
        normalLoopTail.append(sizeRanges['gene'][i])

with open('normalLoopTail.pkl', 'wb') as fp:
    pickle.dump(normalLoopTail, fp)

with open('longLoop.pkl', 'wb') as fp:
    pickle.dump(longLoop, fp)

with open('longTail.pkl', 'wb') as fp:
    pickle.dump(longTail, fp)

with open('longLoopTail.pkl', 'wb') as fp:
    pickle.dump(longLoopTail, fp)

loopCut = np.mean(sizeRanges['icl3']) + np.std(sizeRanges['icl3']) 
tailCut = np.mean(sizeRanges['ctail']) + np.std(sizeRanges['ctail'])

print loopCut, tailCut
