# calculate AT/A/T pct of I bp upstream and downstream of theta
# find TATAWAWR upstream of theta, with 0, 1, 2 mismatches
# generate 200 random samples, each with length 2*I, among all genes with length >700 and do the same statistics



# cd /proj/strahllb/users/Jie/Stephen/ATpct
# module add bedtools
# bsub Rscript ../ATpct_random_sample.R 200 100 121 s

library(Biostrings)

args = commandArgs(trailingOnly = TRUE)
options(scipen = 999)


#args[1] == number of random samples 
#args[2] == I bp surround theta
#args[3] == 28 if for the 28 high genes
        #== 92 if for the 96 genes
        #== 121 if for all 121 genes or any other situation
#args[4] == 's' if with -s option in getfasta
#        == 'us' if without -s
   
result.mat = read.table('/proj/strahllb/users/Jie/Stephen/cryp_result/all_output_rm_badmatch.txt', header = T, sep='\t')

K= as.numeric(args[1]) # of ramdom samples

I = as.numeric(args[2]) # I bps upstream and downstream of theta
 

High = result.mat[result.mat$Group=='High',]

n=nrow(High)



if (as.character(args[3]) == '28') {
print('how many genes:')
print(args[3])
#the gene name with 4+ in all the three time points
high_names = c("YAL043C","YBL035C","YDR097C","YDR108W","YDR181C","YDR325W",
               "YER114C","YGL120C","YGL142C","YGR089W","YHR120W","YHR197W",
               "YJL198W","YJR041C","YKR024C","YLL015W","YLL035W","YLR260W",
               "YMR114C","YMR259C","YNL023C","YNL262W","YOL075C","YOL078W",
               "YOL138C","YOR093C" ,"YOR207C","YPR032W")
kp_idx = c()
for (i in 1:n) {
  if (High$Gene_Name[i] %in% high_names) {
    kp_idx = c(kp_idx, i)
  }
}

High = High[kp_idx,]
n=nrow(High)
}


if (as.character(args[3]) == '92'){
print('how many genes:')
print(args[3])
# # the gene list which in 120min antisense data, the avg read of left of theta < avg right
rm_names= c("YMR149W", "YBR153W", "YNL095C", "YMR114C", "YEL055C", "YNL020C", 
            "YOR093C", "YBL008W", "YPR161C", "YDR181C", "YJL198W", "YDL148C", 
            "YNR011C", "YFL002C", "YJL109C", "YJL129C", "YHR114W", "YDR180W", 
            "YML061C", "YDL112W", "YAL043C", "YPL085W", "YDR325W", "YLR386W", 
            "YLR272C", "YJR041C", "YBR245C", "YNL250W", "YCL014W")
rm_idx = c()
for (i in 1:n) {
  if (High$Gene_Name[i] %in% rm_names) {
    rm_idx = c(rm_idx, i)
  }
}

High = High[-rm_idx,]
n=nrow(High)
}


# making bed files of th surrounding regions --------------------------------------------------------
# Notice! only generate these file when asking for AT percentage, so need to run AT before A or T
th120 = High$START + as.numeric(as.character(High$Best_Theta )) - 1

for (i in 1:n) {
  if (High$STRAND[i] == '+') {
    th120[i] = th120[i] -80
  } else if (High$STRAND[i] == '-') {
    th120[i] = th120[i] +80
  }
}


bed120 = cbind(as.character(High$CHROM), th120-I, th120+I, as.character(High$Gene_Name), rep(1, length(th120)), as.character(High$STRAND))

write.table(bed120, paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes.bed'), quote=F, row.names=F, col.names=F, sep='\t')

if (args[4] == 's') {
  system(paste0('bsub bedtools getfasta -tab -s -fi /proj/dllab/Austin/sacCer3/sacCer3.fa -bed theta_surrounding_', I, 'bp_', nrow(High),'genes.bed -fo theta_surrounding_', I, 'bp_', nrow(High),'genes_Strand.txt'))
} else if (args[4] =='us') {
  system(paste0('bsub bedtools getfasta -tab -fi /proj/dllab/Austin/sacCer3/sacCer3.fa -bed theta_surrounding_', I, 'bp_', nrow(High),'genes.bed -fo theta_surrounding_', I, 'bp_', nrow(High),'genes_noStrand.txt'))
}

pool = which(result.mat$STOP - result.mat$START > 700) 
for (k in 1:K) {
  
  set.seed(k)
  rand_idx = sample(pool, n)
  Rand_gene = result.mat[rand_idx,]
  
  rand_th = rep(NA, n)
  for (i in 1:n){
    set.seed(k)
    rand_th[i] = sample((Rand_gene$START[i]+150) : (Rand_gene$STOP[i]-150), 1) }
  rand_bed = cbind(as.character(Rand_gene$CHROM), rand_th-I, rand_th+I, as.character(Rand_gene$Gene_Name), rep(1, n), as.character(Rand_gene$STRAND))
  write.table(rand_bed, paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.bed'), quote=F, row.names=F, col.names=F, sep='\t')
  if (args[4] == 's') {
    system(paste0('bsub bedtools getfasta -tab -s -fi /proj/dllab/Austin/sacCer3/sacCer3.fa -bed theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.bed -fo theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '_Strand.txt'))
  } else if (args[4] =='us') {
  system(paste0('bsub bedtools getfasta -tab -fi /proj/dllab/Austin/sacCer3/sacCer3.fa -bed theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.bed -fo theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '_noStrand.txt'))
  }
}
  
  

no_use = read.table('/proj/strahllb/users/Jie/Stephen/yeast_0min_senseStrand_full_genes_filledNAs_named.coord/yeast_WT_0min_Rep3_senseStrand_full_genes_filledNAs_named.coord', sep="\t", row.names=1, header=T)
  




ATleft = rep(NA, K+1)
ATright= rep(NA, K+1)

Aleft = rep(NA, K+1)
Aright= rep(NA, K+1)

Tleft = rep(NA, K+1)
Tright= rep(NA, K+1)

mm0 = rep(NA, K+1)
mm1 = rep(NA, K+1)
mm2 = rep(NA, K+1)

# read in bedtools result
if (args[4] == 's') {
 sequ = read.table(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_Strand.txt'), sep='', header=F)
 } else if (args[4] =='us') {
 sequ = read.table(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_noStrand.txt'), sep='', header=F)
 }
sequ = as.character(sequ[,2])



score = rep(NA, n)
mis.match = matrix(NA, n, 6)
colnames(mis.match) = c('0 mismatch', '1 mismatch', '2 mismatch', 'a.s. 0 mismatch', 'a.s 1 mismatch', 'a.s 2 mismatch')
for (i in 1:n){
mis.match0= matchPattern("TATAWAWR", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=0, min.mismatch = 0)
mis.match1= matchPattern("TATAWAWR", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=1, min.mismatch = 1)
mis.match2= matchPattern("TATAWAWR", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=2, min.mismatch = 1)

as.mis.match0= matchPattern("YWTWTATA", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=0, min.mismatch = 0)
as.mis.match1= matchPattern("YWTWTATA", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=1, min.mismatch = 1)
as.mis.match2= matchPattern("YWTWTATA", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=2, min.mismatch = 1)


mis.match[i,] = c(length(mis.match0@ranges), length(mis.match1@ranges), length(mis.match2@ranges), length(as.mis.match0@ranges), length(as.mis.match1@ranges), length(as.mis.match2@ranges))
score[i] = 1*length(mis.match2@ranges) + 5*length(mis.match1@ranges) + 10*length(mis.match0@ranges)
}

mis.match.logic = mis.match >=1

mm0[K+1] = sum(mis.match.logic[,1]) + sum(mis.match.logic[,4])

mm1[K+1] = sum(mis.match.logic[,2]) + sum(mis.match.logic[,5])

mm2[K+1] = sum(mis.match.logic[,3]) + sum(mis.match.logic[,6])








window_size = 10


sequ2 = matrix(NA, n,2*I)
name=c(); 
for (j in I:1) {name = c(name, paste0('upstream',j))}
for (j in 0:(I-1)) {name = c(name, paste0('downstream',j))}
colnames(sequ2)=name



for (i in 1:n) {
  sequ2[i,]=strsplit(sequ[i],'')[[1]]
}


ATpct = rep(NA, 2*I-window_size+1)
Apct = rep(NA, 2*I-window_size+1)
Tpct = rep(NA, 2*I-window_size+1)
Gpct = rep(NA, 2*I-window_size+1)
Cpct = rep(NA, 2*I-window_size+1)


for (j in 1:(2*I-window_size+1)){  
    ATpct[j]= sum(sequ2[,j:(j+window_size-1)] =='A' | sequ2[,j:(j+window_size-1)] =='T')/(n*window_size)
    Apct[j]= sum(sequ2[,j:(j+window_size-1)] =='A')/(n*window_size)
    Tpct[j]= sum(sequ2[,j:(j+window_size-1)] =='T')/(n*window_size)  
    Gpct[j]= sum(sequ2[,j:(j+window_size-1)] =='G')/(n*window_size)
    Cpct[j]= sum(sequ2[,j:(j+window_size-1)] =='C')/(n*window_size) 
}

ATleft[K+1] = mean(ATpct[1:I])
ATright[K+1] = mean(ATpct[(I+1): length(ATpct)])

Aleft[K+1] = mean(Apct[1:I])
Aright[K+1] = mean(Apct[(I+1): length(Apct)])

Tleft[K+1] = mean(Tpct[1:I])
Tright[K+1] = mean(Tpct[(I+1): length(Tpct)])

png(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_ws', window_size, '_', args[4],'.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

plot((-I):(I-window_size+1-1), Apct, type='l', 
     ylim=c(0.1,0.4),  xlab='position to theta', ylab='A and T pct', main=paste0(n, ' genes'))
lines((-I):(I-window_size+1-1), Tpct, col='red')
lines((-I):(I-window_size+1-1), (Tpct+Apct), col='blue')
lines((-I):(I-window_size+1-1), Gpct, col='green4')
lines((-I):(I-window_size+1-1), Cpct, col='orange')
abline(v=0)
legend('topright', c('A pct', 'T pct', 'AT pct', 'G pct', 'C pct'), lty=1, col=c('black', 'red', 'blue', 'green4', 'orange'))
dev.off()

  

for (k in 1:K) {

	if (args[4] == 's') {
	  sequ = read.table(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '_Strand.txt'), sep='', header=F)
	} else if (args[4] =='us') {  
	  sequ = read.table(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '_noStrand.txt'), sep='', header=F)
	}
  sequ = as.character(sequ[,2])
    
	score = rep(NA, n)
	mis.match = matrix(NA, n, 6)
	colnames(mis.match) = c('0 mismatch', '1 mismatch', '2 mismatch', 'a.s. 0 mismatch', 'a.s 1 mismatch', 'a.s 2 mismatch')
	for (i in 1:n){
	mis.match0= matchPattern("TATAWAWR", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=0, min.mismatch = 0)
	mis.match1= matchPattern("TATAWAWR", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=1, min.mismatch = 1)
	mis.match2= matchPattern("TATAWAWR", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=2, min.mismatch = 1)

	as.mis.match0= matchPattern("YWTWTATA", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=0, min.mismatch = 0)
	as.mis.match1= matchPattern("YWTWTATA", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=1, min.mismatch = 1)
	as.mis.match2= matchPattern("YWTWTATA", substr(DNAString(sequ[i]), 1, 200), fixed=FALSE, max.mismatch=2, min.mismatch = 1)


	mis.match[i,] = c(length(mis.match0@ranges), length(mis.match1@ranges), length(mis.match2@ranges), length(as.mis.match0@ranges), length(as.mis.match1@ranges), length(as.mis.match2@ranges))
	score[i] = 1*length(mis.match2@ranges) + 5*length(mis.match1@ranges) + 10*length(mis.match0@ranges)
	}

	mis.match.logic = mis.match >=1

	mm0[k] = sum(mis.match.logic[,1]) + sum(mis.match.logic[,4])

	mm1[k] = sum(mis.match.logic[,2]) + sum(mis.match.logic[,5])

	mm2[k] = sum(mis.match.logic[,3]) + sum(mis.match.logic[,6])


 
  sequ2 = matrix(NA, n, 2*I)
  name=c(); 
  for (j in I:1) {name = c(name, paste0('upstream',j))}
  for (j in 0:(I-1)) {name = c(name, paste0('downstream',j))}
  colnames(sequ2)=name
    
  
  for (i in 1:n) {
    sequ2[i,]=strsplit(sequ[i],'')[[1]]
  }
  
  
  ATpct = rep(NA, 2*I-window_size+1)
  Apct = rep(NA, 2*I-window_size+1)
  Tpct = rep(NA, 2*I-window_size+1)
  Gpct = rep(NA, 2*I-window_size+1)
  Cpct = rep(NA, 2*I-window_size+1)

  
  
  for (j in 1:(2*I-window_size+1)){
    ATpct[j]= sum(sequ2[,j:(j+window_size-1)] =='A' | sequ2[,j:(j+window_size-1)] =='T')/(n*window_size)
    Apct[j]= sum(sequ2[,j:(j+window_size-1)] =='A')/(n*window_size)
    Tpct[j]= sum(sequ2[,j:(j+window_size-1)] =='T')/(n*window_size)
    Gpct[j]= sum(sequ2[,j:(j+window_size-1)] =='G')/(n*window_size)
    Cpct[j]= sum(sequ2[,j:(j+window_size-1)] =='C')/(n*window_size) 
  }
  
  
  ATleft[k] = mean(ATpct[1:I])
  ATright[k] = mean(ATpct[(I+1): length(ATpct)])
  
  Aleft[k] = mean(Apct[1:I])
  Aright[k] = mean(Apct[(I+1): length(Apct)])
  
  Tleft[k] = mean(Tpct[1:I])
  Tright[k] = mean(Tpct[(I+1): length(Tpct)])
  

png(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_ws', window_size,'_seed_',k, '_', args[4], '.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

plot((-I):(I-window_size+1-1), Apct, type='l', 
     ylim=c(0.1,0.4),  xlab='position to theta', ylab='A and T pct', main=paste0(n, ' genes'))
lines((-I):(I-window_size+1-1), Tpct, col='red')
lines((-I):(I-window_size+1-1), (Tpct+Apct), col='blue')
lines((-I):(I-window_size+1-1), Gpct, col='green4')
lines((-I):(I-window_size+1-1), Cpct, col='orange')
abline(v=0)
legend('topright', c('A pct', 'T pct', 'AT pct', 'G pct', 'C pct'), lty=1, col=c('black', 'red', 'blue', 'green4', 'orange'))
dev.off()



}


ATratio = ATleft/ATright
ATprob = mean(ATratio[1:K] < ATratio[K+1])

Aratio = Aleft/Aright
Aprob = mean(Aratio[1:K] < Aratio[K+1])

Tratio = Tleft/Tright
Tprob = mean(Tratio[1:K] < Tratio[K+1])

mm0prob = mean(mm0[1:K] < mm0[K+1])
mm1prob = mean(mm1[1:K] < mm1[K+1])
mm2prob = mean(mm2[1:K] < mm2[K+1])

png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'genes_ATpct_ws', window_size, '_', args[4], '.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(ATratio[1:K], breaks=100, freq = FALSE,
     xlab=paste0( 'AT pct of upstream / downstream'),
     main = paste0('Histogram of AT percentage ratio \n Number of random sample: ', K, '\n p-value = ', 1-ATprob))
lines(density(ATratio[1:K]), lwd = 2)
abline(v=ATratio[K+1], lwd=2, col='red')

dev.off() 

png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'genes_Apct_ws', window_size, '_', args[4],'.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(Aratio[1:K], breaks=100, freq = FALSE,
     xlab=paste0( 'A pct of upstream / downstream'),
     main = paste0('Histogram of A percentage ratio \n Number of random sample: ', K, '\n p-value = ', 1-Aprob))
lines(density(Aratio[1:K]), lwd = 2)
abline(v=Aratio[K+1], lwd=2, col='red')

dev.off()  

png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'genes_Tpct_ws', window_size, '_', args[4],'.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(Tratio[1:K], breaks=100, freq = FALSE,
     xlab=paste0( 'T pct of upstream / downstream'),
     main = paste0('Histogram of T percentage ratio \n Number of random sample: ', K, '\n p-value = ', 1-Tprob))
lines(density(Tratio[1:K]), lwd = 2)
abline(v=Tratio[K+1], lwd=2, col='red')

dev.off() 


png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'genes_TATAbox_0_mismatch', '_', args[4],'.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(mm0[1:K], breaks=100, freq = FALSE,
     xlab=paste0( 'Number of 0 mismatch'),
     main = paste0('Histogram of number of 0 mismatch \n Number of random sample: ', K, '\n p-value = ', 1-mm0prob))
lines(density(mm0[1:K]), lwd = 2)
abline(v=mm0[K+1], lwd=2, col='red')

dev.off()


png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'genes_TATAbox_1_mismatch', '_', args[4],'.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(mm1[1:K], breaks=100, freq = FALSE,
     xlab=paste0( 'Number of 1 mismatch'),
     main = paste0('Histogram of number of 1 mismatch \n Number of random sample: ', K, '\n p-value = ', 1-mm1prob))
lines(density(mm1[1:K]), lwd = 2)
abline(v=mm1[K+1], lwd=2, col='red')

dev.off()


png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'genes_TATAbox_2_mismatch', '_', args[4],'.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(mm2[1:K], breaks=100, freq = FALSE,
     xlab=paste0( 'Number of 2 mismatch'),
     main = paste0('Histogram of number of 2 mismatch \n Number of random sample: ', K, '\n p-value = ', 1-mm2prob))
lines(density(mm2[1:K]), lwd = 2)
abline(v=mm2[K+1], lwd=2, col='red')

dev.off()
  