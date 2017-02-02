cd /proj/strahllb/users/Jie/Stephen/ATpct
module add bedtools
bsub Rscript ../ATpct_random_sample.R 200

args = commandArgs(trailingOnly = TRUE)

#args[1] == number of random samples 


# making bed files of th surrounding regions --------------------------------------------------------

result.mat = read.table('/proj/strahllb/users/Jie/Stephen/cryp_result/all_output_rm_badmatch.txt', header = T, sep='\t')


I = 140 # I bps upstream and downstream of theta
 
K= as.numeric(args[1]) # of ramdom samples

High = result.mat[result.mat$Group=='High',]

n=nrow(High)

# 
# # the gene name with 4+ in all the three time points
# high_names = c("YAL043C","YBL035C","YDR097C","YDR108W","YDR181C","YDR325W",
#                "YER114C","YGL120C","YGL142C","YGR089W","YHR120W","YHR197W",
#                "YJL198W","YJR041C","YKR024C","YLL015W","YLL035W","YLR260W",
#                "YMR114C","YMR259C","YNL023C","YNL262W","YOL075C","YOL078W",
#                "YOL138C","YOR093C" ,"YOR207C","YPR032W")
# kp_idx = c()
# for (i in 1:n) {
#   if (High$Gene_Name[i] %in% high_names) {
#     kp_idx = c(kp_idx, i)
#   }
# }
# 
# High = High[kp_idx,]
# n=nrow(High)
# 
# # the gene list which in 120min antisense data, the avg read of left of theta < avg right
# rm_names= c("YMR149W", "YBR153W", "YNL095C", "YMR114C", "YEL055C", "YNL020C", 
#             "YOR093C", "YBL008W", "YPR161C", "YDR181C", "YJL198W", "YDL148C", 
#             "YNR011C", "YFL002C", "YJL109C", "YJL129C", "YHR114W", "YDR180W", 
#             "YML061C", "YDL112W", "YAL043C", "YPL085W", "YDR325W", "YLR386W", 
#             "YLR272C", "YJR041C", "YBR245C", "YNL250W", "YCL014W")
# rm_idx = c()
# for (i in 1:n) {
#   if (High$Gene_Name[i] %in% rm_names) {
#     rm_idx = c(rm_idx, i)
#   }
# }
# 
# High = High[-rm_idx,]
# n=nrow(High)

th120 = High$START + as.numeric(as.character(High$Best_Theta )) - 1

bed120 = cbind(as.character(High$CHROM), th120-I, th120+I, as.character(High$Gene_Name), rep(1, length(th120)), as.character(High$STRAND))

write.table(bed120, paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes.bed'), quote=F, row.names=F, col.names=F, sep='\t')

system(paste0('bsub bedtools getfasta -s -tab -fi /proj/dllab/Austin/sacCer3/sacCer3.fa -bed theta_surrounding_', I, 'bp_', nrow(High),'genes.bed -fo theta_surrounding_', I, 'bp_', nrow(High),'genes.txt'))


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
  system(paste0('bsub bedtools getfasta -s -tab -fi /proj/dllab/Austin/sacCer3/sacCer3.fa -bed theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.bed -fo theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.txt'))
  
}
  
  

random = read.table('/proj/strahllb/users/Jie/Stephen/yeast_0min_senseStrand_full_genes_filledNAs_named.coord/yeast_WT_0min_Rep3_senseStrand_full_genes_filledNAs_named.coord', sep="\t", row.names=1, header=T)
  
  



left = rep(NA, K+1)
right= rep(NA, K+1)

# read in bedtools result
sequ = read.table(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes.txt'), sep='', header=F)
sequ = as.character(sequ[,2])


window_size = 1


sequ2 = matrix(NA, n,2*I)
name=c(); 
for (j in I:1) {name = c(name, paste0('upstream',j))}
for (j in 0:(I-1)) {name = c(name, paste0('downstream',j))}
colnames(sequ2)=name



for (i in 1:n) {
  sequ2[i,]=strsplit(sequ[i],'')[[1]]
}


ATpct = rep(NA, 2*I-window_size+1)


for (j in 1:(2*I-window_size+1)){
  ATpct[j]= sum(sequ2[,j:(j+window_size-1)] =='A' | sequ2[,j:(j+window_size-1)] =='T')/(n*window_size)
}

left[K+1] = mean(ATpct[1:I])
right[K+1] = mean(ATpct[(I+1): length(ATpct)])

png(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

plot((-I):(I-window_size+1-1), ATpct, type='l', 
     ylim=c(0.4,0.8),  xlab='position to theta', ylab='AT pct', main=paste0(n, ' genes'))
abline(h=0.6)
# abline(v=100)
# abline(v=-100)
abline(v=0)

dev.off()

  

for (k in 1:K) {
  sequ = read.table(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.txt'), sep='', header=F)
  sequ = as.character(sequ[,2])
  
  
  window_size = 1
  
  
  sequ2 = matrix(NA, n, 2*I)
  name=c(); 
  for (j in I:1) {name = c(name, paste0('upstream',j))}
  for (j in 0:(I-1)) {name = c(name, paste0('downstream',j))}
  colnames(sequ2)=name
  
  
  
  for (i in 1:n) {
    sequ2[i,]=strsplit(sequ[i],'')[[1]]
  }
  
  
  ATpct = rep(NA, 2*I-window_size+1)
  
  
  for (j in 1:(2*I-window_size+1)){
    ATpct[j]= sum(sequ2[,j:(j+window_size-1)] =='A' | sequ2[,j:(j+window_size-1)] =='T')/(n*window_size)
  }
  
  
  left[k] = mean(ATpct[1:I])
  right[k] = mean(ATpct[(I+1): length(ATpct)])
  
#   png(paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k, '.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  
#   
#   plot((-I):(I-window_size+1-1), ATpct, type='l', 
#        ylim=c(0.4,0.8),  xlab='position to theta', ylab='AT pct', main=paste0('theta_surrounding_', I, 'bp_', nrow(High),'genes_seed_', k))
#   abline(h=0.6)
#   # abline(v=100)
#   # abline(v=-100)
#   abline(v=0)
#   
#   dev.off()
}


ratio = left/right
prob = mean(ratio[1:K] < ratio[K+1])

print('left average ATpct:')
left
print('righ average ATpct:')
right

png(paste0('hist_theta_surrounding_', I, 'bp_', nrow(High),'_genes.png'), width = 8*300,  height = 5*300,  res = 300, pointsize = 8)  

hist(ratio[1:K], breaks=100, freq = FALSE,
     xlab='ATpct of upstream / downstream',
     main = paste0('Histogram of AT percentage ratio \n Number of random sample: ', K, '\n p-value = ', 1-prob))
lines(density(ratio[1:K]), lwd = 2)
abline(v=ratio[K+1], lwd=2, col='red')

dev.off() 


print('the probability of AT ratio is higher in the upstream is:')
print(prob)

  