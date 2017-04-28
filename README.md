# PhosCoFinder
A python based tool, paired with http://tools.vai.org/phoscofinder/, to rapidly identify potential phosphorylation codes for arrestin recruitment. The accompanying manuscript, Structural Identification of Phosphorylation Codes ofr Arrestin Recruitment by G protein-Coupled Receptors (Zhou et al., 2017), can be found at TBD.

## Scanning Sequences
PhosCoFinder.py is launched with the following arguments:
```
python phoscofinder.py {search name} {sequence location}
```
An example quicksearch of all known annotated GPCRs from uniprot can be launched with:
```
python phoscofinder.py uniprot-GPCR-ctail fasta/uniprot/uniprot-gpcrs.fasta
```

## Supplied Datasets
PhosCoFinder is provided with 3 sources of pre-formated data to scan: 
1. [GPCRdb Class divided C-Tail sequences](/fasta/gpcrdb)
2. [UniProt](/fasta/uniprot)
  * All annotated GPCR C-termianl tail sequences
  * All annotated GPCR ICL3 sequences
3. [Known bArr1/2 interactors from the Proteomics Standard Initiative Common QUery InterfaCe (PSICQUIC)](/fasta/psicquic) 

Additional source links and instructions on how to format `phoscofinder` sequence files can be found within each databases respected folder.
