`uniprot.py` is provided as a conversion tool from single/multi protein .txt files obtained from [Uniprot](http://http://www.uniprot.org/) to PhosCoFinder compatible sequence files. Initial serach strings are provided below

* [Verified human GPCRs search string](http://www.uniprot.org/uniprot/?query=gpcr&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score)
* [Verified, human non-GPCR membrane proteins obtained from uniprot search string](http://www.uniprot.org/uniprot/?query=annotation%3A%28type%3Atopo_dom+Cytoplasmic%29+AND+annotation%3A%28type%3Atransmem%29+NOT+go%3A%22g+protein%22+AND+reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+NOT+keyword%3A%22G-protein+coupled+receptor+%5BKW-0297%5D%22&sort=score)

To convert sequences, `uniprot.py` requires the following arguments:
```
python uniprot.py {GPCR/non} {ICL3/cTail} uniprot.txt > uniprot.fasta
```

To extract ICL3 and C-tail sequences:
```
python uniprot.py GPCR ICL3 uniprot-gpcrs.txt > uniprot-ICL3-gpcrs.txt.fasta
python uniprot.py GPCR cTail uniprot-gpcrs.txt > uniprot-cTail-gpcrs.txt.fasta
```

Non-GPCR membrane proteins can be extracted similary:
```
python uniprot.py non cTail uniprot-non-gpcrs.txt > uniprot-non-gpcrs.txt.fasta
```

