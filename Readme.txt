performing test run:


bsub Rscript ~/determine_cryptic_jumps_and_rates_faster.R ~/single_gene_seq.txt ~/output.txt ~/loss.txt 0

# args[1] == inputfile, args[2] == outfile1, arg[3] == outfile2,
# args[4] == 0

make sure the following 3 R packages are installed: 
lpSolve
Matrix
methods


format of single_gene_seq.txt file:
1st line: WT group gene info
2nd line: SET2 deleted group gene info
1st column: Gene name
2nd column: Chromosom number
3rd column: CDS start position
4th column: CDS end position
5th column: 1
6th column: strand infomation
7th column and after: reads starting from CDS position 1 to the end of CDS (No nan allowed)

format of the output.txt file:
a list of numbers with the following information:
genename, gene length, cryptic initiation site, y, z, MSRLD.

formate of the loss.txt file:
a list of number with the following information:
MSRL value of each position of the gene CDS (except for the starting and ending 150bp respectively).

rule of thumb: 
High cryptic initiation: MSRLD >=4
Intermediate cryptic initiation: 2 <= MSRLD < 4
Low cryptic initiation: 0 <= MSRLD < 2
Not quite trustworthy (bad fitting): Min(MSRL) > 15
