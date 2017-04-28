Formated sequence files obtained from the [GPCRdb](http://gpcrdb.org/) only contail the C-tail region and thus are numbered incorrectly.

To generate PhosCoFinder compatible files, simply download fasta files from the GPCRdb and convert using `convert.sh`

Ex:
```
bash convert.sh classA.fasta > classA.fasta.oneLine
```
