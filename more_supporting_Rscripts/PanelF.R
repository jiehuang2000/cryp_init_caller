out = read.table('~/all_output_rm_badmatch_yz.txt', header=T, sep='\t')

data = read.table('~/Downloads/yeast_SET2del_120min_Avg_antisenseStrand_full_genes_filledNAs_named.coord', header=T, sep='\t',row.names=1)

data = cbind(rownames(data), data)
colnames(data)[1] = 'Gene_Name'
data = merge(out, data, by='Gene_Name')


gene_start = which(colnames(data) == "Position1")
L=ncol(data)

high = data[data$Group == 'High'|data$Group == 'Medium',]


as.avg=rep(NA,nrow(high))
for (i in 1:nrow(high)) {
  as = high[i, gene_start:L]
  as = as[!is.na(as)]
  th = high$Best_Theta[i] 
  as.avg[i] = mean(as[1:th])
}

as.avg.hm = as.avg

#as.avg[as.avg==0] =1e-3


# lm = lm(y ~ x + high$Best_z_120 + high$Best_Theta_120)
# summary(lm)
# abline(lm$coefficients)
# 
# library(MASS)
# bla <-boxcox(lm,plotit=T,lambda=seq(0.1, 0.5, by=0.1))
# 
# 
# cor(as.avg, high$Best_y_120)
# cor(as.avg, high$Best_z_120)
# cor(high$Best_y_120, high$Best_z_120)
# cor(as.avg, high$Loss_Diff_120)
# cor(high$Best_z_120, high$Loss_Diff_120)
# cor(high$Best_y_120, high$Loss_Diff_120)


 



WT = read.table('~/Downloads/yeast_WT_120min_Avg_antisenseStrand_full_genes_filledNAs_named.coord', header=T, sep='\t',row.names=1)

WT = cbind(rownames(WT), WT)
colnames(WT)[1] = 'Gene_Name'
WT = merge(out, WT, by='Gene_Name')

high.WT = WT[WT$Group == 'High'|WT$Group == 'Medium',]

as.avg.WT=rep(NA,nrow(high.WT))
for (i in 1:nrow(high.WT)) {
  as = high.WT[i, gene_start:L]
  as = as[!is.na(as)]
  th = high.WT$Best_Theta[i] - 50
  as.avg.WT[i] = mean(as[1:th])
}
# 
# temp=as.avg - as.avg.WT
# hist(temp)
# cor(temp, high$Best_y_120)
# cor(temp, high$Best_z_120)
# cor(temp, high$Loss_Diff_120)
# 
# plot(high$Loss_Diff_120, temp)
# lm=lm(temp~high$Loss_Diff_120)
# summary(lm)
# abline(lm$coefficients)
# bla <-boxcox(lm,plotit=T,lambda=seq(0.1, 0.5, by=0.1))




x=log2(high$Best_y_120)
#y=as.avg-as.avg.WT
y=log2((as.avg+0.1)/(as.avg.WT+0.1))
sum(y>0 & x<0)/length(as.avg)

pdf('panelF.pdf', width=5, height=5)
par(mar=c(4,6,4,4))
smoothScatter(x, y, 
              xlab='Predicted full length mRNA expression ratio\n(SET2del/WT)',
              ylab='Antisense expression ratio\n(SET2del/WT)',
              colramp = colorRampPalette(c("white", 'cyan3',  'lightgreen', 'yellow','orange','red')),transformation = function(x) (x^0.9)) 
abline(v=0,h=0)
lcb5=which(high$Gene_Name=='YLR260W')
avo1=which(high$Gene_Name=='YOL078W')
points(x[lcb5],y[lcb5], pch=4, col='red')
text(x[lcb5]  ,y[lcb5]+1.5, 'lcb5')

points(x[avo1],y[avo1], pch=4, col='red')
text(x[avo1]  ,y[lcb5]+1.5, 'avo1')
legend('topleft', legend=paste0(round(sum(y>0 & x<0)/length(as.avg),4)*100, '%'))

# 
dev.off() 
