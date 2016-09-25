the determine_cryptic_jumps_and_rates.R algorithm test run using the given single_gene_seq.txt

bsub Rscript ~/determine_cryptic_jumps_and_rates_faster.R ~/single_gene_seq.txt ~/output.txt ~/loss.txt 0

# args[1] == inputfile, args[2] == outfile1:genename, gene.len, theta, y, z, Diff , arg[3] == outfile, loss value of all sites
# args[4] == 0/1, if we want stepsize=1, args[4] = 0, if want stepsize=100, args[4] == 1


