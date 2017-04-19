##############################################
# Zhou et al. Cell 2016 -- Parker de Waal 2017
# Simple conversion of FASTA files downloaded from GPCR database (gpcrdb.org)
#
# Usage ex: bash process.sh classA.fasta classA.fasta.oneline
##############################################

sed 's/\-//g' $1 > tmp
awk '/>/{gene=$1; getline; if(length($1) > 0){print gene, " 1",length($1),$1}}' tmp > $2
sed -i 's/>//' $2
rm tmp
