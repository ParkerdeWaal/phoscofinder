##############################################
# Zhou et al. Cell 2016 -- Parker de Waal 2017
# Simple conversion of FASTA files downloaded from GPCR database (gpcrdb.org)
#
# Usage ex: bash process.sh classA.fasta classA.fasta.oneline
##############################################

awk '/>/{gene=$1; getline; print gene, $1}' $1 > $2
