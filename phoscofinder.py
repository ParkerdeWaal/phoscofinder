#!/usr/bin/env python

###################################################
# Find Regex matches for phosphotail motifs PxPxxP and PxxPxxP from GPCRdb fasta files
# Zhou Et al. Cell (2017) -- by Parker de Waal
# 
# Usage ex: python phoscofinder.py A fasta/classA.fasta.online
#                               (class)  
#
#
#
##################################################

import sys,os
import re
import numpy as np
import numpy.lib.recfunctions as rfn
from time import gmtime, strftime
from collections import Counter

# define class and fasta file
inClass = sys.argv[1]
inFile = sys.argv[2]
sequences = open(inFile, "r")
folder = "./"+strftime("%Y-%m-%d-%H-%M-%S-", gmtime()) + str(inClass) + "/"
os.mkdir(folder)

# Regex rules
# 1. PxPxxP or PxxPxxP where P1 or P2 are either S/T. P3 can be STE or D
# 2. Last 4 amino acids cannot contain a proline

regexShort = re.compile(r'(S|T)(?=\w{1}[S|T][^\W_P]{2}[S|T|E|D])')
regexLong = re.compile(r'(S|T)(?=\w{2}[S|T][^\W_P]{2}[S|T|E|D])')

# define partial short motifs
ps1 = re.compile(r'(S|T)(?=\w{1}[S|T][^\W_P]{2}[^\W_STEDP])')
ps2 = re.compile(r'(S|T)(?=\w{1}[^\W_STP][^\W_P]{2}[S|T|E|D])')
ps3 = re.compile(r'([^\W_ST])(?=\w{1}[S|T][^\W_P]{2}[S|T|E|D])')

# define partial long motifs
pl1 = re.compile(r'(S|T)(?=\w{2}[S|T][^\W_P]{2}[^\W_STEDP])')
pl2 = re.compile(r'(S|T)(?=\w{2}[^\W_STP][^\W_P]{2}[S|T|E|D])')
pl3= re.compile(r'[^\W_ST](?=\w{2}(S|T)[^\W_P]{2}[S|T|E|D])')


# Define array type and initialize array
geneMotif = np.dtype([('gene', 'a11'),('start', 'i'),('stop', 'i'),('match','a7')])
geneOnly = np.dtype([('gene', 'a11')])
countArr = np.dtype([('gene', 'a11'),('count', 'i')])

foundShort = np.empty([0,2],dtype=geneMotif)
foundLong = np.empty([0,2],dtype=geneMotif)
foundNone = np.empty([0,2],dtype=geneOnly)
foundPartial = np.empty([0,2],dtype=geneMotif)

totalReceptors = 0

# Loop over sequences and identify motifs
for line in sequences:
	totalReceptors +=1
	foundMotif = 0
	geneName = re.sub(r'\>', '',line.split()[0])
	sequence = line.split()[3][10:] #ignore first 10 resdiues as they MAY be in H8
	print geneName
	if int(line.split()[1]) != 1:
		offset = int(line.split()[1])
	else:
		offset = 1 + 1
	for match in re.finditer(regexShort,sequence):
		partialSeq = sequence[(match.start()):(match.end()+5)]
		tmpStart = int(match.start())+offset
		tmpEnd = match.end()+offset+4
		foundShort = np.append(foundShort, np.array([(geneName,match.start()+offset,match.end()+offset+4, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"S",tmpStart,tmpEnd,partialSeq)
		foundMotif +=1

	for match in re.finditer(regexLong,sequence):
		partialSeq = sequence[(match.start()):(match.end()+6)]
		foundLong = np.append(foundLong, np.array([(geneName,match.start()+offset,match.end()+offset+5, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"L",match.start()+offset,match.end()+offset+5,partialSeq)
		foundMotif +=1

	for match in re.finditer(ps1,sequence):
		partialSeq = sequence[(match.start()):(match.end()+5)]
		foundPartial = np.append(foundPartial, np.array([(geneName,match.start()+offset,match.end()+offset+4, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"PS",match.start()+offset,match.end()+offset+4,partialSeq)
		foundMotif +=1

	for match in re.finditer(ps2,sequence):
		partialSeq = sequence[(match.start()):(match.end()+5)]
		foundPartial = np.append(foundPartial, np.array([(geneName,match.start()+offset,match.end()+offset+4, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"PS",match.start()+offset,match.end()+offset+4,partialSeq)
		foundMotif +=1

	for match in re.finditer(ps3,sequence):
		partialSeq = sequence[(match.start()):(match.end()+5)]
		foundPartial = np.append(foundPartial, np.array([(geneName,match.start()+offset,match.end()+offset+4,partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"PS",match.start()+offset,match.end()+offset+4,partialSeq)
		foundMotif +=1

	for match in re.finditer(pl1,sequence):
		partialSeq = sequence[(match.start()):(match.end()+6)]
		foundPartial = np.append(foundPartial, np.array([(geneName,match.start()+offset,match.end()+offset+5, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"PL",match.start()+offset,match.end()+offset+5,partialSeq)
		foundMotif +=1

	for match in re.finditer(pl2,sequence):
		partialSeq = sequence[(match.start()):(match.end()+6)]
		foundPartial = np.append(foundPartial, np.array([(geneName,match.start()+offset,match.end()+offset+5, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"PL",match.start()+offset,match.end()+offset+5,partialSeq)
		foundMotif +=1

	for match in re.finditer(pl3,sequence):
		partialSeq = sequence[(match.start()):(match.end()+6)]
		foundPartial = np.append(foundPartial, np.array([(geneName,match.start()+offset,match.end()+offset+5, partialSeq)], dtype=geneMotif))
		print "%11s %1s %4d %4d %17s" % (geneName,"PL",match.start()+offset,match.end()+offset+5,partialSeq)
		foundMotif +=1

	if (foundMotif == 0):
		geneName = re.sub(r'\>', '',line.split()[0])
		foundNone = np.append(foundNone, np.array([(geneName)], dtype=geneOnly))

numbShortUnique = np.count_nonzero(np.unique(foundShort['gene'])) # number of short receptors
numbLongUnique = np.count_nonzero(np.unique(foundLong['gene'])) # number of long receptors
numbInterUnique = np.count_nonzero(np.unique(np.intersect1d(foundShort['gene'], foundLong['gene']))) # total number of receptors containing complete code
numbWith = np.count_nonzero(np.unique(np.union1d(foundShort['gene'], foundLong['gene'])))
numbWithout = np.count_nonzero(np.unique(foundNone['gene']))
numbPartial = totalReceptors - numbWith - numbWithout


print "Short:%4d   Long:%4d   Intercept:%4d   With:%4d    Partial:%4d    Without:%4d    total:%4d" %(numbShortUnique,numbLongUnique,numbInterUnique,numbWith,numbPartial,numbWithout,totalReceptors)


interData = np.unique(np.intersect1d(foundShort['gene'], foundLong['gene']))


shortCount = Counter(foundShort['gene'])
longCount = Counter(foundLong['gene'])
partialCount = Counter(foundPartial['gene'])

shortCount = np.array(sorted(shortCount.items(), key=lambda pair: pair[1], reverse=True), dtype=countArr)
longCount = np.array(sorted(longCount.items(), key=lambda pair: pair[1], reverse=True),dtype=countArr)
partialCount = np.array(sorted(partialCount.items(), key=lambda pair: pair[1], reverse=True),dtype=countArr)

joined = rfn.join_by('gene', shortCount, longCount, jointype='outer',usemask=False)
joined = rfn.join_by('gene', joined, partialCount, jointype='outer',usemask=False)
joined = rfn.join_by('gene', joined, foundNone, jointype='outer',usemask=False)


for x in range(0,joined.size):
	if joined['count1'][x] == 999999:
		joined['count1'][x] = 0
        if joined['count2'][x] == 999999:
                joined['count2'][x] = 0
        if joined['count'][x] == 999999:
                joined['count'][x] = 0

fileData = folder + inClass + ".codes.csv"
fileNameFound = folder + inClass + ".found.csv"
fileNameNone = folder + inClass + ".none.csv"
 

np.savetxt(fileData,joined, delimiter=',',fmt=('%11s ',' %5i',' %5i',' %5i'), header='uniprotID, shortCodes,longCodes, partialCodes,',comments='')
fo =  np.unique(np.append(foundShort['gene'],foundLong['gene']))
fo =  np.unique(np.append(foundPartial['gene'],fo))
np.savetxt(fileNameFound,fo,delimiter=',',fmt=('%11s'), header='uniprotID')
np.savetxt(fileNameNone,foundNone,delimiter=',',fmt=('%11s'), header='uniprotID')
