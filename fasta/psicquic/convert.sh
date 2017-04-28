## Working script to convert psique database files to phoscofinder sequences
rm arrb12.txt tmp ids unique


cat arrb1.txt arrb2.txt > arrb12.txt

touch tmp
awk '{print $1}' arrb12.txt >> tmp
awk '{print $2}' arrb12.txt >> tmp

touch ids
awk -F':' '{print $2}' tmp >> ids

touch unique
sort ids | uniq >> unique

python downloader.py unique > uniprot.txt

python uniprot.py non cTail uniprot.txt > knownInteractors.fasta


sort knownInteractors.fasta | uniq > tmp; mv tmp knownInteractors.fasta
