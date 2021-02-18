# Geneious Prime Steps

1. Create a folder for your project or dataset for eg. Haiti.
2. Create a sub folder called Individual or Pooled based on your sample type.
3. Import all bam files here.
4. Geneious will automatically identify which sample contains reads for which genes.
5. Using Geneious workflow, separate these into 4 different folders, one for each gene (dhfr, mdr1, dhps, crt)
6. Using workflow again, for dhfr and mdr1 "map to reference" gene sequences and select save consensus sequence; for dhps and crt "map to reference" the CDS sequence of the respective genes and select save consensus sequence.
7. In each gene folder, select all the newly assembled sequences, click of "Text View" and them "Save as *.txt".
8. This will generate a single document containing all the sequences assembled to the selected gene.
9. Use this as one of the inputs to the calc_coverage.py which will process this document and use each sample sequence to obtain start and end positions for coverage calculations.

![](https://github.com/MGupta313/Codon-coverage/blob/master/Geneious/Folders.png)

![picture alt](http://via.placeholder.com/200x150 "Title is optional")