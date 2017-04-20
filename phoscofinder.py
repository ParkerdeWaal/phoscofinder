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
	sequence = line.split()[3]
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


fileNameShort = folder + inClass + ".short.csv"
fileNameLong = folder + inClass + ".long.csv"
fileNameInter = folder + inClass + ".inter.csv"
fileNameNone = folder + inClass + ".none.csv"
fileNamePartial = folder + inClass + ".partial.csv"

np.savetxt(fileNameShort,foundShort, delimiter=',',fmt=('%11s ',' %5i',' %5i', ' %16s'), header='GeneID,Sequence',comments='')
np.savetxt(fileNameLong,foundLong, delimiter=',',fmt=('%11s ',' %5i',' %5i', ' %16s'), header='GeneID,Sequence',comments='')
np.savetxt(fileNamePartial,foundPartial, delimiter=',',fmt=('%11s ',' %5i',' %5i', ' %16s'), header='GeneID,Sequence',comments='')
np.savetxt(fileNameInter,interData, delimiter=',',fmt=('%11s'), header='GeneID',comments='')
np.savetxt(fileNameNone,foundNone, delimiter=',',fmt=('%11s'), header='GeneID',comments='')

shortCount = Counter(foundShort['gene'])
longCount = Counter(foundLong['gene'])
shortlongCount = Counter(foundShort['gene']) + Counter(foundLong['gene'])
partialCount = Counter(foundPartial['gene'])

shortCount = np.array(sorted(shortCount.items(), key=lambda pair: pair[1], reverse=True), dtype=countArr)
longCount = np.array(sorted(longCount.items(), key=lambda pair: pair[1], reverse=True),dtype=countArr)
shortlongCount = np.array(sorted(shortlongCount.items(), key=lambda pair: pair[1], reverse=True),dtype=countArr)
partialCount = np.array(sorted(partialCount.items(), key=lambda pair: pair[1], reverse=True),dtype=countArr)

np.savetxt(folder + inClass + ".shortCount.csv",shortCount, delimiter=',',fmt=('%11s','%5i'), header='GeneID,count',comments='')
np.savetxt(folder + inClass + ".longCount.csv",longCount, delimiter=',',fmt=('%11s','%5i'), header='GeneID,count',comments='')
np.savetxt(folder + inClass + ".shortlongCount.csv",shortlongCount, delimiter=',',fmt=('%11s','%5i'), header='GeneID,count',comments='')
np.savetxt(folder + inClass + ".partialCount.csv",partialCount, delimiter=',',fmt=('%11s','%5i'), header='GeneID,count',comments='')
