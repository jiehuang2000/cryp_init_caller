
# import cryp init data ---------------------------------------------------
gene_info = read.table('~/Dropbox/JieHuang/plot_raw_files/all_6692_genes_info_name_chrom_start_strop_strand.txt', header=T, sep='\t')


d60min = read.table('~/Dropbox/JieHuang/plot_raw_files/final_matrix_name_of_all_genes_60min_.out', header=T, sep='\t')
head(d60min)
colnames(d60min)=c('Gene_Name', 'Note_60', 'Best_Theta_60', 'Best_y_60', 'Best_z_60', 'Loss_Diff_60', 'Loss_min_60')


d30min = read.table('~/Dropbox/JieHuang/plot_raw_files/final_matrix_name_of_all_genes_30min_.out', header=T, sep='\t')
head(d30min)
colnames(d30min)=c('Gene_Name', 'Note_30', 'Best_Theta_30', 'Best_y_30', 'Best_z_30', 'Loss_Diff_30', 'Loss_min_30')



d0min = read.table('~/Dropbox/JieHuang/plot_raw_files/final_matrix_name_of_all_genes_0min_.out', header=T, sep='\t')
####test data
#d0min = read.table('~/Downloads/final_matrix_name_of_all_genes_0min_WT1_mut2_.out', header=T, sep='\t')
head(d0min)
colnames(d0min)=c('Gene_Name', 'Note_0', 'Best_Theta_0', 'Best_y_0', 'Best_z_0', 'Loss_Diff_0', 'Loss_min_0')

d120min = read.table('~/Dropbox/JieHuang/plot_raw_files/final_matrix_name_of_all_genes_120min_.out', header=T, sep='\t')
head(d120min)
colnames(d120min)=c('Gene_Name', 'Note_120', 'Best_Theta_120', 'Best_y_120', 'Best_z_120', 'Loss_Diff_120', 'Loss_min_120')

dim(d60min)
dim(d30min)
dim(d0min)
dim(d120min)



summary(d60min$Note_60)[1:5]
summary(d30min$Note_30)[1:5]
summary(d0min$Note_0)[1:5]
summary(d120min$Note_120)[1:5]




# plot loss histogram -----------------------------------------------------

# png(paste0('Histogram of MMSRL for each gene.png'), width = 5*600,  height = 5*600,  res = 600, pointsize = 8)  
pdf(paste0('Histogram of MMSRL for each gene.pdf'), width = 5,  height = 5)

hist(d0min$Loss_min_0, breaks=seq(0,90), col=rgb(0,0,0,0.1), ylim=c(0,0.3),freq = FALSE , 
     main='Histogram of MMSRL for each gene',
     xlab='MMSRL', ylab='Density')
lines(density(d0min$Loss_min_0, na.rm=T), lwd = 2)

hist(d30min$Loss_min_30, breaks=seq(0,30), col=rgb(1,0,0,0.3) , add=T,freq = FALSE)
lines(density(d30min$Loss_min_30, na.rm=T), lwd = 2, col=rgb(1,0,0))

hist(d60min$Loss_min_60, breaks=seq(0,30), col=rgb(0,1,0,0.3) , add=T,freq = FALSE)
lines(density(d60min$Loss_min_60, na.rm=T), lwd = 2, col=rgb(0,1,0))

hist(d120min$Loss_min_120, breaks=seq(0,30), col=rgb(0,0,1,0.2) , add=T,freq = FALSE)
lines(density(d120min$Loss_min_120, na.rm=T), lwd = 2, col=rgb(0,0,1))

cutoff = mean(c(d30min$Loss_min_30, d60min$Loss_min_60, d120min$Loss_min_120), na.rm=T) + 1.96* sd(c(d30min$Loss_min_30, d60min$Loss_min_60, d120min$Loss_min_120), na.rm=T)

abline(v=cutoff, lty=3, col='orange')


legend('topright', c('0min','30min', '60min', '120min', 'cutoff'), col = c('black', 'red', 'green', 'blue', 'orange'), lty = c(1,1,1,1,3))
dev.off()

sum(d60min$Loss_min_60 >= cutoff, na.rm=T)
sum(d30min$Loss_min_30 >= cutoff, na.rm=T)
sum(d0min$Loss_min_0 >= cutoff, na.rm=T)
sum(d120min$Loss_min_120 >= cutoff, na.rm=T)



# remove bad match & add category -----------------------------------------

######60min################
temp=matrix(NA, nrow(d60min),1)
for (i in 1:nrow(d60min)){
  #print(i)
  if (is.na(d60min$Loss_min_60[i]) | is.na(d60min$Loss_Diff_60[i])) {
    temp[i,1] = as.character(d60min$Note_60[i])
    d60min[i, 3:7] = rep(NA, 5)
  }else if (d60min$Loss_min_60[i] >= cutoff) {
    temp[i,1] = 'bad_match'
    d60min[i, 3:7] = rep(NA, 5)
  }else if (d60min$Loss_Diff_60[i] >= 4 ) {
    temp[i,1] = 'High'
  }else if (d60min$Loss_Diff_60[i] >= 2 & d60min$Loss_Diff_60[i] < 4 ) {
    temp[i,1] = 'Medium'
  }else {
    temp[i,1] = 'Low'
  }
}
Category_60 = as.factor(temp)
d60min = cbind(d60min, Category_60)

write.table(d60min, file='final_matrix_name_of_all_genes_60min_bad_match_as_NA.out', quote=F, sep="\t", col.names=T, row.names=F)

######30min################
temp=matrix(NA, nrow(d30min),1)
for (i in 1:nrow(d30min)){
  #print(i)
  if (is.na(d30min$Loss_min_30[i]) | is.na(d30min$Loss_Diff_30[i])) {
    temp[i,1] = as.character(d30min$Note_30[i])
    d30min[i, 3:7] = rep(NA, 5)
  }else if (d30min$Loss_min_30[i] >= cutoff) {
    temp[i,1] = 'bad_match'
    d30min[i, 3:7] = rep(NA, 5)
  }else if (d30min$Loss_Diff_30[i] >= 4 ) {
    temp[i,1] = 'High'
  }else if (d30min$Loss_Diff_30[i] >= 2 & d30min$Loss_Diff_30[i] < 4 ) {
    temp[i,1] = 'Medium'
  }else {
    temp[i,1] = 'Low'
  }
}
Category_30 = as.factor(temp)
d30min = cbind(d30min, Category_30)
write.table(d30min, file='final_matrix_name_of_all_genes_30min_bad_match_as_NA.out', quote=F, sep="\t", col.names=T, row.names=F)

######120min################
temp=matrix(NA, nrow(d120min),1)
for (i in 1:nrow(d120min)){
  #print(i)
  if (is.na(d120min$Loss_min_120[i]) | is.na(d120min$Loss_Diff_120[i])) {
    temp[i,1] = as.character(d120min$Note_120[i])
    d120min[i, 3:7] = rep(NA, 5)
  }else if (d120min$Loss_min_120[i] >= cutoff) {
    temp[i,1] = 'bad_match'
    d120min[i, 3:7] = rep(NA, 5)
  }else if (d120min$Loss_Diff_120[i] >= 4 ) {
    temp[i,1] = 'High'
  }else if (d120min$Loss_Diff_120[i] >= 2 & d120min$Loss_Diff_120[i] < 4 ) {
    temp[i,1] = 'Medium'
  }else {
    temp[i,1] = 'Low'
  }
}
Category_120 = as.factor(temp)
d120min = cbind(d120min, Category_120)
write.table(d120min, file='final_matrix_name_of_all_genes_120min_bad_match_as_NA.out', quote=F, sep="\t", col.names=T, row.names=F)

######0min################
temp=matrix(NA, nrow(d0min),1)
for (i in 1:nrow(d0min)){
  #print(i)
  if (is.na(d0min$Loss_min_0[i]) | is.na(d0min$Loss_Diff_0[i])) {
    temp[i,1] = as.character(d0min$Note_0[i])
    d0min[i, 3:7] = rep(NA, 5)
  }else if (d0min$Loss_min_0[i] >= cutoff) {
    temp[i,1] ='bad_match'
    d0min[i, 3:7] = rep(NA, 5)
  }else if (d0min$Loss_Diff_0[i] >= 4 ) {
    temp[i,1] = 'High'
  }else if (d0min$Loss_Diff_0[i] >= 2 & d0min$Loss_Diff_0[i] < 4 ) {
    temp[i,1] = 'Medium'
  }else {
    temp[i,1] = 'Low'
  }
}
Category_0 = as.factor(temp)
d0min = cbind(d0min, Category_0)
write.table(d0min, file='final_matrix_name_of_all_genes_0min_bad_match_as_NA.out', quote=F, sep="\t", col.names=T, row.names=F)
#write.table(d0min, file='final_matrix_name_of_all_genes_0min_WT2_mut2_bad_match_as_NA.out', quote=F, sep="\t", col.names=T, row.names=F)


dim(d60min)
dim(d30min)
dim(d0min)
dim(d120min)

summary(Category_0)
summary(Category_30)
summary(Category_60)
summary(Category_120)

sum(!is.na(d0min$Loss_Diff_0))
sum(!is.na(d30min$Loss_Diff_30))
sum(!is.na(d60min$Loss_Diff_60))
sum(!is.na(d120min$Loss_Diff_120))



# remove NA genes ---------------------------------------------------------
# rmID=c()
# for (i in 1:nrow(d60min)) {
#   if (is.na(d60min$Loss_Diff_60[i])) {
#     rmID=c(rmID,i)
#   }
# }
# d60min=d60min[-rmID,]
# 
# 
# rmID=c()
# for (i in 1:nrow(d30min)) {
#   if (is.na(d30min$Loss_Diff_30[i])) {
#     rmID=c(rmID,i)
#   }
# }
# d30min=d30min[-rmID,]
# 
# 
# rmID=c()
# for (i in 1:nrow(d0min)) {
#   if (is.na(d0min$Loss_Diff_0[i])) {
#     rmID=c(rmID,i)
#   }
# }
# d0min=d0min[-rmID,]
# 
# 
# rmID=c()
# for (i in 1:nrow(d120min)) {
#   if (is.na(d120min$Loss_Diff_120[i])) {
#     rmID=c(rmID,i)
#   }
# }
# d120min=d120min[-rmID,]
# 
# dim(d60min)
# dim(d30min)
# dim(d0min)
# dim(d120min)



# flipped results ---------------------------------------------------------
d60min_flipped = read.table('~/Dropbox/JieHuang/plot_raw_files/flipped_final_matrix_name_of_all_genes_60min.out', header=T, sep='\t')
head(d60min_flipped)
colnames(d60min_flipped)=c('Gene_Name', 'Gene_Length', 'Best_Theta_60', 'Best_y_60', 'Best_z_60', 'Loss_Diff_60')

d30min_flipped = read.table('~/Dropbox/JieHuang/plot_raw_files/flipped_final_matrix_name_of_all_genes_30min.out', header=T, sep='\t')
head(d30min_flipped)
colnames(d30min_flipped)=c('Gene_Name', 'Gene_Length', 'Best_Theta_30', 'Best_y_30', 'Best_z_30', 'Loss_Diff_30')

d0min_flipped = read.table('~/Dropbox/JieHuang/plot_raw_files/flipped_final_matrix_name_of_all_genes_0min.out', header=T, sep='\t')
head(d0min_flipped)
colnames(d0min_flipped)=c('Gene_Name', 'Gene_Length', 'Best_Theta_0', 'Best_y_0', 'Best_z_0', 'Loss_Diff_0')

d120min_flipped = read.table('~/Dropbox/JieHuang/plot_raw_files/flipped_final_matrix_name_of_all_genes_120min.out', header=T, sep='\t')
head(d120min_flipped)
colnames(d120min_flipped)=c('Gene_Name', 'Gene_Length', 'Best_Theta_120', 'Best_y_120', 'Best_z_120', 'Loss_Diff_120')

dim(d60min_flipped)
dim(d30min_flipped)
dim(d0min_flipped)
dim(d120min_flipped)





# merge data --------------------------------------------------------------



d= merge(d0min,d30min, by='Gene_Name')
dim(d)
d=merge(d,d60min, by='Gene_Name')
dim(d)
d=merge(d,d120min, by='Gene_Name')
dim(d)

d= merge(gene_info, d, by.x='COORDID', by.y='Gene_Name')







# categorize genes into groups: High Medium Low Bad ------------------------
# Bad: if any of the 3 time points (30min, 60min, 120min) is NA
# High: if any of the 3 time points is >=4
# Low: if all of the 3 time points is <2
# Mediuma: anything else 


loss_diff = cbind(d$Loss_Diff_0, d$Loss_Diff_30, d$Loss_Diff_60, d$Loss_Diff_120)
rownames(loss_diff)=as.character(d[,1])
head(loss_diff)
colnames(loss_diff) = c('0min', '30min', '60min', '120min')
# make loss_diff = 0 to loss_diff = 1e-9 for the later on log2 transformation
loss_diff[!is.na(loss_diff) & loss_diff==0] = 1e-9


Group = rep(NA, nrow(d))
for (i in 1:nrow(d)){
  if ( sum(is.na(loss_diff[i,2:4])) >= 1 ) {
    Group[i] = 'Bad'
  #} else if (loss_diff[i,2] >= 4 & loss_diff[i,3] >= 4 & loss_diff[i,4] >= 4 ){
  } else if (loss_diff[i,2] >= 4 | loss_diff[i,3] >= 4 | loss_diff[i,4] >= 4 ){
    Group[i] = 'High'
  } else if (loss_diff[i,2] < 2 & loss_diff[i,3] < 2 & loss_diff[i,4] < 2) {
    Group[i] = 'Low'
  } else {Group[i] = 'Medium'}
}

Group=as.factor(Group)
summary(Group)

#d= cbind(Group, d)




# remove data with NA for any 30 60 120min time point and make heatmap ---------------
# 
# NAid=c()
# for (i in 1:nrow(loss_diff)){
#   if ( sum(is.na(loss_diff[i,2:4])) >= 1 ) {NAid=c(NAid,i)}
#   for (j in 1:4){
#     if (is.na(loss_diff[i,j])) {loss_diff[i,j] = -2}
#   }
# }
# 
# loss_diff_rmNA=loss_diff[-NAid,]
# dim(loss_diff_rmNA)
# colnames(loss_diff_rmNA) = c('0min', '30min', '60min', '120min')
# 
# 
# all0id=c()
# for (i in 1:nrow(loss_diff)){
#   for (j in 1:4){
#     if (is.na(loss_diff[i,j])) {loss_diff[i,j] = -2}  
#   }
#   if (sum(loss_diff[i,] < 1) == 4) {all0id=c(all0id,i)}
# }
# 
# loss_diff_rmNA_lessthan1_read=loss_diff[-all0id,]
# dim(loss_diff_rmNA_lessthan1_read)
# colnames(loss_diff_rmNA_lessthan1_read) = c('0min', '30min', '60min', '120min')



# all0id=c()
# for (i in 1:nrow(loss_diff)){
#   if ( sum(is.na(loss_diff[i,])) == 4 ) {all0id=c(all0id,i)}
#   for (j in 1:4){
#     if (is.na(loss_diff[i,j])) { loss_diff[i,j] = -2
#     }else if (loss_diff[i,j] >=4) { loss_diff[i,j] = 4
#     }else if (loss_diff[i,j] >=2 ) { loss_diff[i,j] = 2
#     }else {loss_diff[i,j] = 0}                                 
#   }
#   if (sum(loss_diff[i,] == 0)==4) {all0id=c(all0id,i)}
# }
# 
# loss_diff=loss_diff[-all0id,]
# dim(loss_diff)
# colnames(loss_diff) = c('0min', '30min', '60min', '120min')



# Make heatmap ------------------------------------------------------------

library(gplots)
options(expressions=500000)


loss_diff_no0 = loss_diff[,2:4]
colnames(loss_diff_no0) = c('30min', '60min', '120min')




my_palette <- colorRampPalette(c("white", "blue"))(n = 40)
#png(paste0('cryp_result_rm_bad_match_exclude_0min_ge4_log2_transformed.png'), width = 5*600,  height = 8*600,  res = 600, pointsize = 8)  
pdf(paste0('cryp_result_rm_bad_match_exclude_0min_ge4_log2_transformed.pdf'), width=6, height=10) 
heatmap.2(data.matrix(log2(loss_diff_no0[Group == 'High',])),dendrogram="row", 
#heatmap.2(data.matrix(log2(loss_diff[Group == 'High',])),dendrogram="row",          
          col=my_palette, 
          breaks=seq(0,4,length.out=41),
          trace = 'none',
          Colv=FALSE , 
          srtCol = 360, cexCol=2, cexRow = .5, adjCol = c(0.5,1),
          #labRow = NA,
          density.info = "none" , keysize=1, 
          lmat=rbind(c(4, 3), c(2, 1)), lhei=c(1, 8), lwid=c(1,4),
          key.xlab='Log2 transformed \nMSRL difference',
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 1)),
          main=paste0('High Chance \n' ,nrow(loss_diff_no0[Group == 'High', ]), ' genes '))

dev.off()

my_palette <- c(colorRampPalette(c("white", "yellow"))(n = 10), 
                colorRampPalette(c("yellow", "skyblue1"))(n = 11)[2:11],
                colorRampPalette(c("skyblue1", "blue"))(n = 21)[2:21])
#png(paste0('cryp_result_rm_bad_match_exclude_0min_ge2_log2_transformed.png'), width = 5*600,  height = 8*600,  res = 600, pointsize = 8)  
pdf(paste0('cryp_result_rm_bad_match_exclude_0min_ge2_log2_transformed.pdf'), width=6, height=10) 
heatmap.2(data.matrix(log2(loss_diff_no0[Group == 'High'| Group == 'Medium',])),dendrogram="row", 
          col=my_palette, breaks=seq(0,4,length.out=41),
          trace = 'none',
          Colv=FALSE , 
          srtCol = 360, cexCol=2, cexRow = .5, adjCol = c(0.5,1),
          labRow = NA,
          density.info = "none" , keysize=1, 
          lmat=rbind(c(4, 3), c(2, 1)), lhei=c(1, 8), lwid=c(1,4),
          key.xlab='Log2 transformed \nMSRL difference',
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 1)),
          main=paste0('High & Medium Chance \n' ,nrow(loss_diff_no0[Group == 'High' | Group == 'Medium', ]), ' genes '), cex=.8)

dev.off()
 
#png(paste0('cryp_result_rm_bad_match_exclude_0min_ge0_log2_transformed.png'), width = 5*600,  height = 8*600,  res = 600, pointsize = 8) 
pdf(paste0('cryp_result_rm_bad_match_exclude_0min_ge0_log2_transformed.pdf'), width=6, height=10) 
heatmap.2(data.matrix(log2(loss_diff_no0[Group == 'High'| Group == 'Medium' | Group == 'Low',])),dendrogram="row", 
          col=my_palette, breaks=seq(0,4,length.out=41),
          trace = 'none',
          Colv=FALSE , 
          srtCol = 360, cexCol=2, cexRow = .5, adjCol = c(0.5,1),
          labRow = NA,
          density.info = "none" , keysize=1, 
          lmat=rbind(c(4, 3), c(2, 1)), lhei=c(1, 6), lwid=c(1,4),
          key.xlab='Log2 transformed \nMSRL difference',
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 1)),
          main=paste0('All Genes \n' ,nrow(loss_diff_no0[Group == 'High' | Group == 'Medium'| Group == 'Low', ]), ' genes '))

dev.off()

# heatmap.2(data.matrix(log2(loss_diff_no0[Group == 'High',])),dendrogram="row", 
#           col=my_palette, breaks=seq(0,4,length.out=41),
#           Colv=FALSE, density.info = "none" , keysize=0.8, 
#           srtCol = 360, cexCol=2, cexRow = .5, adjCol = c(0.5,1),
#           key.par=list(mgp=c(1.5, 0.5, 0),
#                        mar=c(2.5, 2.5, 1, 0)),
#           main=paste0(nrow(loss_diff_no0[Group == 'High', ]), ' genes'))
# 
# dev.off()

# #http://stackoverflow.com/questions/5687891/r-how-do-i-display-clustered-matrix-heatmap-similar-color-patterns-are-grouped
# k = kmeans(loss_diff, 5)
# dfc <- cbind(loss_diff, id=seq(nrow(loss_diff)), cluster=k$cluster)
# dfc$idsort <- dfc$id[order(dfc$cluster)]
# dfc$idsort <- order(dfc$idsort)
# dfm <- melt(dfc, id.vars=c("id", "idsort"))
# ggplot(dfm, aes(x=variable, y=idsort)) + geom_tile(aes(fill=value))
# 
 


# make best_th.txt file ---------------------------------------------------


best_th = cbind(d$Best_Theta_0, d$Best_Theta_30, d$Best_Theta_60, d$Best_Theta_120)
rownames(best_th)=as.character(d[,1])
head(best_th)

#loss_diff = round(d[,c(6,11,16,21)])
loss_diff = cbind(d$Loss_Diff_0, d$Loss_Diff_30, d$Loss_Diff_60, d$Loss_Diff_120)
rownames(loss_diff)=as.character(d[,1])
head(loss_diff)
colnames(loss_diff) = c('0min', '30min', '60min', '120min')
# make loss_diff = 0 to loss_diff = 1e-9 for the later on log2 transformation
loss_diff[!is.na(loss_diff) & loss_diff==0] = 1e-9


#find best the for all time points. 
th=c()

for (i in 1:nrow(best_th)){
  if (sum(is.na(loss_diff[i,2:4])) ==3) {
    th=c(th, NA)
    
  }else{
    idx= which(loss_diff[i,2:4] ==max(loss_diff[i,2:4], na.rm=T))[1]
    th=c(th, best_th[i,idx+1])
  }
}

# th=round(apply(best_th, 1, mean, na.rm=T))
th = cbind(as.character(d$COORDID), th, as.character(Group))
colnames(th) = c('Gene_Name', 'Best_Theta', "Group")
write.table(th, file='best_th_Group.txt', sep="\t", col.names=T, row.names=F, quote=F)

#output.matrix = d[, -c(7,9,10,12, 14,16,17,19, 21,23,24,26, 28,30,31,33)]
output.matrix = d
output.matrix=merge(th, output.matrix , by.x='Gene_Name', by.y='COORDID')
write.table(output.matrix, file='all_output_rm_badmatch.txt', ,sep='\t', quote=F, row.names=F, col.names=T)







# plot of 'Jump point position variations for 3 time points.pdf' ----------

thetas = output.matrix[output.matrix$Group=='High', c(which(colnames(output.matrix)=='Best_Theta_30'),
                                             which(colnames(output.matrix)=='Best_Theta_60'),
                                             which(colnames(output.matrix)=='Best_Theta_120'))]

thetas2 = output.matrix[output.matrix$Group=='Medium', c(which(colnames(output.matrix)=='Best_Theta_30'),
                                                      which(colnames(output.matrix)=='Best_Theta_60'),
                                                      which(colnames(output.matrix)=='Best_Theta_120'))]
png('Jump point position variations for 3 time points.png', width = 5*600,  height = 8*600,  res = 600, pointsize = 8)

#pdf('Jump point position variations for 3 time points.pdf', width=6, height=8)
par(mfrow=c(2,1))
plot(1:nrow(thetas),thetas[,1], col='blue', cex=0.7,
     ylim=c(0, max(thetas)+100), ylab='Jump point position',
     xlab='121 genes in group "High"',
     main='Jump point position variations for 3 time points')
points(1:nrow(thetas),thetas[,2],pch=3,col='orange', cex=0.7)
points(1:nrow(thetas),thetas[,3],pch=4,col='green4', cex=0.7)

for (i in 1:nrow(thetas)){
lines(c(i,i), c(max(thetas[i,]),min(thetas[i,])),  col='yellow', lty=3)
}
legend('topright', c('30min', '60min', '120min'), col=c('blue', 'orange', 'green4'), pch=c(1,3,4), cex=0.7)

plot(1:nrow(thetas2),thetas2[,1], col='blue', cex=0.7,
     ylim=c(0, max(thetas2)+100), ylab='Jump point position',
     xlab='318 genes in group "Medium"',
     main='Jump point position variations for 3 time points')
points(1:nrow(thetas2),thetas2[,2],pch=3,col='orange', cex=0.7)
points(1:nrow(thetas2),thetas2[,3],pch=4,col='green4', cex=0.7)

for (i in 1:nrow(thetas2)){
  lines(c(i,i), c(max(thetas2[i,]),min(thetas2[i,])),  col='yellow', lty=3)
}
legend('topright', c('30min', '60min', '120min'), col=c('blue', 'orange', 'green4'), pch=c(1,3,4), cex=0.7)
dev.off()

thetas_sd = apply(thetas,1,sd)

par(new=T)
plot(1:length(thetas_sd), thetas_sd, axes=F, type='line')
abline(h=200)
sum(thetas_sd <=200)/length(thetas_sd)

png(paste0('Sorted delta plot of four time points.png'), width = 5*600,  height = 5*600,  res = 600, pointsize = 8)  


d0minL = sort(d0min$Loss_Diff_0,decreasing = TRUE)
d30minL = sort(d30min$Loss_Diff_30,decreasing = TRUE)
d60minL = sort(d60min$Loss_Diff_60,decreasing = TRUE)
d120minL = sort(d120min$Loss_Diff_120,decreasing = TRUE)

d0minL_flipped = sort(d0min_flipped$Loss_Diff_0, decreasing = TRUE)
d30minL_flipped = sort(d30min_flipped$Loss_Diff_30, decreasing = TRUE)
d60minL_flipped = sort(d60min_flipped$Loss_Diff_60, decreasing = TRUE)
d120minL_flipped = sort(d120min_flipped$Loss_Diff_120, decreasing = TRUE)


plot(d0minL,type='l',ylim = c(0,15), xlim=c(0,4000),ylab='max loss difference (delta)')
title('Sorted delta plot of four time points')
lines(d30minL,col='red')
lines(d60minL,col='green')
lines(d120minL,col='blue')

lines(d0minL_flipped, lty=3)
lines(d30minL_flipped, lty=3, col='red')
lines(d60minL_flipped, lty=3, col='green')
lines(d120minL_flipped, lty=3, col='blue')

abline(h=4, lty=2, col='orange')
abline(h=2, lty=2, col='orange')
legend('topright', c('0min','30min', '60min', '120min','flipped'), col = c('black', 'red', 'green', 'blue','grey'), lty = c(1,1,1,1,3))

dev.off()


require(lattice)
require(ggplot2)

#png(paste0('pairwise comparison of Delta value for four time points.png'), width = 5*600,  height = 5*600,  res = 600, pointsize = 8)  
png(paste0('pairwise comparison of Delta value for four time points_0min_WT1_mut2_rm_badmatch.png'), width = 5*600,  height = 5*600,  res = 600, pointsize = 8)  
colnames(loss_diff) = c('0min', '30min', '60min', '120min')
pairs(loss_diff, pch = '.', main = 'pairwise comparison of Delta value for four time points')
dev.off()



pdf('Histogram of length of all genes.pdf',width=6, height=4)
len=output.matrix$STOP - output.matrix$START +1
hist(len, breaks=seq(0,15000,by=200), 
     col=rgb(1,0,0,0.4), border=F,
     xlab='Gene length',
     main='Histogram of length of all genes\n median=1086, mean=1369')

for (i in seq(0,15000,by=1000)) {
  for (j in seq(0, 900, by=100)){
    abline(v=i, lty=3, col='gray96')
    abline(h=j, lty=3, col='gray96')
  }
}
dev.off()


d60min[which(d60min$Gene_Name=="YNL160W"),]
d60min[which(d60min$Gene_Name=="YLR427W"),]
d60min[which(d60min$Gene_Name=="YPR161C"),]
d60min[which(d60min$Gene_Name=="YKL105C"),]

d60min[c(5166,5028,5863,5407),]
# Gene_Name Gene_Length Best_Theta_60 Best_y_60 Best_z_60 Loss_Diff_60 Loss_min_60
# 5166   YNL160W        1067           416 1.3265880 0.2962131     4.327318   15.789550
# 5028   YLR427W        2015           641 0.6975645 0.6893394     7.209892    6.707955
# 5863   YPR161C        1976          1037 0.9380828 0.4216616     5.593152    8.055248
# 5407   YKL105C        3401          2188 1.0055490 0.5893506     4.674433    8.180494



# Influence of gene length_on cryptic initiation --------------------------
len[which(len>6000)] = 6001
max_len = max(len)
output.matrix2 = cbind(len, output.matrix)
window_size=500
#png('Influence of gene length_on cryptic initiation.png', width = 6*600,  height = 4*600,  res = 600, pointsize = 8)  
pdf('Influence of gene length_on cryptic initiation.pdf',width=6, height=4)
par(mar = c(5,5,2,5))
hist(output.matrix2$len[output.matrix2$Group!='Bad'], 
     col=rgb(0,0,0,0.5),
     breaks=seq(700,max_len+window_size,by=window_size),
     xlim=c(700,max_len+window_size),
     include.lowest = TRUE,
     xaxt = "n", 
     border=F, 
     #ylim=c(0,150),
     main='Influence of gene length on cryptic initiation',
     xlab=paste0('Gene Length (window size =', window_size, ')'),
     ylab='Number of genes in each gene length range')
axis(side = 1,at = seq(700,max_len+window_size, by = window_size))
hist(output.matrix2$len[output.matrix2$Group=='High'|output.matrix2$Group=='Medium'], 
     col=rgb(1,0,0,0.7),
     breaks=seq(700,max_len+window_size,by=window_size),
     include.lowest = TRUE,
     xaxt = "n", 
     border=F,
     add =T)

hist(output.matrix2$len[output.matrix2$Group=='High'], 
     col=rgb(0,0,1,1), 
     breaks=seq(700,max_len+window_size,by=window_size),
     include.lowest = TRUE,
     border=F, 
     xaxt = "n",
     add=T)
#rug(jitter(output.matrix2$len[output.matrix2$Group=='High']), col='blue')

ge4=table(cut(output.matrix2$len[output.matrix2$Group=='High'],seq(700,max_len+window_size,by=window_size),include.lowest = TRUE))

ge2 = table(cut(output.matrix2$len[output.matrix2$Group=='High'|output.matrix2$Group=='Medium'],seq(700,max_len+window_size,by=window_size),include.lowest = TRUE))

ge0 =  table(cut(output.matrix2$len[output.matrix2$Group!='Bad'],seq(700,max_len+window_size,by=window_size),include.lowest = TRUE))

par(new=T)
plot(seq(700+window_size/2,max_len,by=window_size), (ge4/ge0)[1:length(seq(700+window_size/2,max_len,by=window_size))], pch=3, axes=F,xlab=NA, ylab=NA, cex=0.5, type='o',
     col=rgb(0,0,1,1), ylim=c(0,0.6), lwd=1, lty=1, xlim=c(700,max_len+window_size))

par(new=T)
plot(seq(700+window_size/2,max_len,by=window_size), (ge2/ge0)[1:length(seq(700+window_size/2,max_len,by=window_size))], pch=4, axes=F, xlab=NA, ylab=NA, cex=0.5, type='o',
     col=rgb(1,0,0,1), ylim=c(0,0.6), lwd=1, lty=1, xlim=c(700,max_len+window_size))

axis(side = 4)
mtext(side = 4, line = 3, 'Proportion')
# legend('topright', 
#        c('High frequency','Medium frequency','Low frequency', 'High proportion','Medium proportion'), 
#        col = c(rgb(0,0,1,0.7),rgb(0,0,0,0.7),rgb(1,0,0,0.7),rgb(0,0,1,1),rgb(0,0,0,1)),
#        pch = c(15, 15, 15, 3, 4), )

legend('topright', 
       c('High ','Medium ','Low ', 'Frequency','Proportion'), 
       col = c(rgb(0,0,1,0.7),rgb(1,0,0,0.7),rgb(0,0,0,1),rgb(0,0,0,1),rgb(0,0,0,1)),
       pch = c(20, 20, 20, 15, 8) )
dev.off()


pdf('Influence of gene length_on cryptic initiation_2.pdf', width = 6,  height = 4)  
#png('Influence of gene length_on cryptic initiation_2.png', width = 6*600,  height = 4*600,  res = 600, pointsize = 8)  
th = as.numeric(as.character(output.matrix2$Best_Theta))
hist((th[output.matrix2$Group=='High'|output.matrix2$Group=='Medium']-150) / (output.matrix2$len[output.matrix2$Group=='High'|output.matrix2$Group=='Medium'] -300), breaks=seq(0,1, by=.05), border=F,
     col=rgb(1,0,0,0.5),xlab='jump point position/gene length', main='Relative Jump point position within a gene')
hist((th[output.matrix2$Group=='High']-150) / (output.matrix2$len[output.matrix2$Group=='High'] -300), breaks=seq(0,1, by=.05), border=F,
     col=rgb(0,0,1,1), add=T)

legend('topright', 
       c('High ','Medium '), 
       col = c(rgb(0,0,1,1),rgb(1,0,0,0.7)),
       pch = c(15, 15), )

dev.off()

