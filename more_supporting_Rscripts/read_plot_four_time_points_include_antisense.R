#testing testing                         
#cd /proj/strahllb/users/Jie/Stephen/cryp_result/
#bsub Rscript ../read_plot_four_time_points_include_antisense.R 10

args = commandArgs(trailingOnly = TRUE)

i = args[1]


jpeg(paste0('gene_', i, '.jpg'), width = 600, height = 1000, units = "px")
par(mfrow=c(4,1))

if (!file.exists(paste0('../yeast_0min_senseStrand_full_genes_filledNAs_named.coord/WT1_mut2/gene_',i, '_WTavg_MUTavg.txt')) || !file.exists(paste0('../yeast_30min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg.txt')) ||!file.exists(paste0('../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg.txt'))||!file.exists(paste0('../yeast_120min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg.txt')) ) {
  stop(paste0('file of gene ', i, 'does not exist!')) 
}


d = read.table(paste0('../yeast_0min_senseStrand_full_genes_filledNAs_named.coord/WT1_mut2/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
antisense.d = read.table(paste0('../yeast_0min_antisenseStrand_full_genes_filledNAs_named.coord/single_gene_sequence/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
info = read.table(paste0('../yeast_0min_senseStrand_full_genes_filledNAs_named.coord/WT1_mut2/gene_',i, '_WTavg_MUTavg_out.txt'), header=F, sep='\t')
loss = read.table(paste0('../yeast_0min_senseStrand_full_genes_filledNAs_named.coord/WT1_mut2/gene_',i, '_WTavg_MUTavg_loss.txt'), header=F, sep='\t')
dim(d)
T=length(d)
max = max(loss)
min = min(loss)
d0 = d[1,7:T]
d1 = d[2,7:T]
as.d0 = antisense.d[1,7:T]
as.d1 = antisense.d[2,7:T]
plot(ksmooth(1:length(d0), d0, "normal",bandwidth = 1), type='l', ylim=c(0,300),ylab='read', xlab='position', main=paste0(toString(info[[1]]), ', gene length: ', info[2],'\n', '0min, loss diff=', round(info[6],2), ' ,change point= ', info[3], ', y=', round(info[4],2), ', z=', round(info[5],2), '\n', 'lossMax=', round(max,2), ', lossMin= ', round(min,2), ', (lossMax-lossMin)/lossMin=', round(info[6]/ min,2)))
lines(ksmooth(1:length(d1), d1, "normal",bandwidth = 1), type='l', col=2)
polygon(c(1:length(as.d0),length(as.d0),1), c(as.d0,0,0), col=rgb(0,0,0, 0.2), border=F)
polygon(c(1:length(as.d1),length(as.d1),1), c(as.d1,0,0), col=rgb(1,0,0, 0.2), border=F)

if (info[6] >=4) {
    color = 'lightskyblue'
} else if (info[6] >=2) {
    color = 'black'
} else if (info[6] >=0) {
    color = 'red'
} else {
    color = 'yellow'
}

abline(v=info[3], col=color)
legend('topleft', c('sense: WT', "sense: SET2 deleted", 'antisense: WT', "antisense: SET2 deleted"), col=c('black', 'red',rgb(0,0,0, 0.2), rgb(1,0,0, 0.2)), pch=c(19,19,19,19))

d = read.table(paste0('../yeast_30min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
antisense.d = read.table(paste0('../yeast_30min_antisenseStrand_full_genes_filledNAs_named.coord/single_gene_sequence/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
info = read.table(paste0('../yeast_30min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg_out.txt'), header=F, sep='\t')
loss = read.table(paste0('../yeast_30min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg_loss.txt'), header=F, sep='\t')
dim(d)
T=length(d)
max = max(loss)
min = min(loss)
d0 = d[1,7:T]
d1 = d[2,7:T]
as.d0 = antisense.d[1,7:T]
as.d1 = antisense.d[2,7:T]
plot(ksmooth(1:length(d0), d0, "normal",bandwidth = 1), type='l', ylim=c(0,300),ylab='read', xlab='position', main=paste0( '30min, loss diff=', round(info[6],2), ' ,change point= ', info[3], ', y=', round(info[4],2), ', z=', round(info[5],2), '\n', 'lossMax=', round(max,2), ', lossMin= ', round(min,2), ', (lossMax-lossMin)/lossMin=', round(info[6]/ min,2)))
lines(ksmooth(1:length(d1), d1, "normal",bandwidth = 1), type='l', col=2)
polygon(c(1:length(as.d0),length(as.d0),1), c(as.d0,0,0), col=rgb(0,0,0, 0.2), border=F)
polygon(c(1:length(as.d1),length(as.d1),1), c(as.d1,0,0), col=rgb(1,0,0, 0.2), border=F)

if (info[6] >=4) {
    color = 'lightskyblue'
} else if (info[6] >=2) {
    color = 'black'
} else if (info[6] >=0) {
    color = 'red'
} else {
    color = 'yellow'
}


abline(v=info[3], col=color)
legend('topleft', c('sense: WT', "sense: SET2 deleted", 'antisense: WT', "antisense: SET2 deleted"), col=c('black', 'red',rgb(0,0,0, 0.2), rgb(1,0,0, 0.2)), pch=c(19,19,19,19))

d = read.table(paste0('../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
antisense.d = read.table(paste0('../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/single_gene_sequence/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
info = read.table(paste0('../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg_out.txt'), header=F, sep='\t')
loss = read.table(paste0('../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg_loss.txt'), header=F, sep='\t')
dim(d)
T=length(d)
max = max(loss)
min = min(loss)
d0 = d[1,7:T]
d1 = d[2,7:T]
as.d0 = antisense.d[1,7:T]
as.d1 = antisense.d[2,7:T]
plot(ksmooth(1:length(d0), d0, "normal",bandwidth = 1), type='l', ylim=c(0,300),ylab='read', xlab='position', main=paste0( '60min, loss diff=', round(info[6],2), ' ,change point= ', info[3], ', y=', round(info[4],2), ', z=', round(info[5],2), '\n', 'lossMax=', round(max,2), ', lossMin= ', round(min,2), ', (lossMax-lossMin)/lossMin=', round(info[6]/ min,2)))
lines(ksmooth(1:length(d1), d1, "normal",bandwidth = 1), type='l', col=2)
polygon(c(1:length(as.d0),length(as.d0),1), c(as.d0,0,0), col=rgb(0,0,0, 0.2), border=F)
polygon(c(1:length(as.d1),length(as.d1),1), c(as.d1,0,0), col=rgb(1,0,0, 0.2), border=F)

if (info[6] >=4) {
    color = 'lightskyblue'
} else if (info[6] >=2) {
    color = 'black'
} else if (info[6] >=0) {
    color = 'red'
} else {
    color = 'yellow'
}


abline(v=info[3], col=color)
legend('topleft', c('sense: WT', "sense: SET2 deleted", 'antisense: WT', "antisense: SET2 deleted"), col=c('black', 'red',rgb(0,0,0, 0.2), rgb(1,0,0, 0.2)), pch=c(19,19,19,19))


d = read.table(paste0('../yeast_120min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
antisense.d = read.table(paste0('../yeast_120min_antisenseStrand_full_genes_filledNAs_named.coord/single_gene_sequence/gene_',i, '_WTavg_MUTavg.txt'), header=F, sep='\t')
info = read.table(paste0('../yeast_120min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg_out.txt'), header=F, sep='\t')
loss = read.table(paste0('../yeast_120min_senseStrand_full_genes_filledNAs_named.coord/cryp_init_result/gene_',i, '_WTavg_MUTavg_loss.txt'), header=F, sep='\t')
dim(d)
T=length(d)
max = max(loss)
min = min(loss)
d0 = d[1,7:T]
d1 = d[2,7:T]
as.d0 = antisense.d[1,7:T]
as.d1 = antisense.d[2,7:T]
plot(ksmooth(1:length(d0), d0, "normal",bandwidth = 1), type='l', ylim=c(0,300),ylab='read', xlab='position', main=paste0( '120min, loss diff=', round(info[6],2), ' ,change point= ', info[3], ', y=', round(info[4],2), ', z=', round(info[5],2), '\n', 'lossMax=', round(max,2), ', lossMin= ', round(min,2), ', (lossMax-lossMin)/lossMin=', round(info[6]/ min,2)))
lines(ksmooth(1:length(d1), d1, "normal",bandwidth = 1), type='l', col=2)
polygon(c(1:length(as.d0),length(as.d0),1), c(as.d0,0,0), col=rgb(0,0,0, 0.2), border=F)
polygon(c(1:length(as.d1),length(as.d1),1), c(as.d1,0,0), col=rgb(1,0,0, 0.2), border=F)

if (info[6] >=4) {
    color = 'lightskyblue'
} else if (info[6] >=2) {
    color = 'black'
} else if (info[6] >=0) {
    color = 'red'
} else {
    color = 'yellow'
}


abline(v=info[3], col=color)

legend('topleft', c('sense: WT', "sense: SET2 deleted", 'antisense: WT', "antisense: SET2 deleted"), col=c('black', 'red',rgb(0,0,0, 0.2), rgb(1,0,0, 0.2)), pch=c(19,19,19,19))


dev.off()


