data= read.table('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_E_heatmap/fix_th_yeast_SET2del_120min_Avg_antisenseStrand_full_genes_filledNAs_named.coord_ge4.txt', sep='\t', header=T) 
dim(data)

th = which(colnames(data) == 'downstream0')
gene_start= which(colnames(data) == 'upstream15000')

n = nrow(data)


idx = c()
idx_rm = c()
for (i in 1:n) {
  left = data[i, gene_start:(th-1)]
  left = left[!is.na(left)]
  right = data[i, th: ncol(data)]
  right = right[!is.na(right)]
  if (mean(left) > mean(right)) {
    idx = c(idx, i)
  } else {
    idx_rm=c(idx_rm,i)
  }
}

maxL = length(left)

data1 = data[idx,]
dim(data1)

RHS = 2000
#data2 = cbind(data1[,(th-maxL-100):(th-1)], matrix(NA, nrow(data1), 30), data1[,(th):(th+RHS)])
data2 = cbind(data[,(th-maxL-100):(th-1)], matrix(NA, nrow(data), 30), data[,(th):(th+RHS)])


library(gplots)
options(expressions=500000)

my_palette <- colorRampPalette(c("white", "black"))(n = 100)
#pdf('High Level (Average left greater than right).pdf',width=8, height=4)
pdf('High Level (Average left greater than right) log2.pdf',width=8, height=4)
#png('High Level (Average left greater than right).png',width=8*600, height=4*600, res=600, pointsize=8)
#png('High Level (Average left greater than right) log2.png',width=8*600, height=4*600, res=600, pointsize=8)

#heatmap.2(data.matrix(data2),
heatmap.2(log2(data.matrix(data2)),
          #breaks=seq(0,100,length.out=101), col=my_palette, na.color="wheat",
          breaks=seq(0,7,length.out=101), col=my_palette, na.color="wheat",
          dendrogram="none", Rowv=FALSE, Colv=FALSE,
          trace="none", 
          #density.info = "none" ,
          labCol=NA,
          cexRow = .3, 
          #key.xlab='Read Depth', 
          key.xlab=paste0(expression(log[2]), '(Read Depth)'),
          main = paste0('High Level (Average left > right) \n', 
                        'Upstream ',maxL+100, ' to Downstream ', RHS, '\n' ,
                        nrow(data2), ' genes')
          )
dev.off()




data3 = cbind(data[,(th-maxL-100):(th-1)], matrix(NA, nrow(data), 30), data[,(th):(th+RHS)])


