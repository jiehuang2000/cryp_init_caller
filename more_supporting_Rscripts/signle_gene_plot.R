
#in terminal

  # cd /proj/strahllb/users/Jie/Stephen/test/
  # mkdir 2834
  # rm 2834/*.txt
  # sed -n 2834p ../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/yeast_WT_60min_Rep1_senseStrand_full_genes_filledNAs_named.coord > 2834/2834_sense_WT1.txt
  # sed -n 2834p ../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/yeast_WT_60min_Rep2_senseStrand_full_genes_filledNAs_named.coord > 2834/2834_sense_WT2.txt
  # sed -n 2834p ../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/yeast_WT_60min_Rep3_senseStrand_full_genes_filledNAs_named.coord > 2834/2834_sense_WT3.txt
  # sed -n 2834p ../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/yeast_SET2del_60min_Rep1_senseStrand_full_genes_filledNAs_named.coord > 2834/2834_sense_MUT1.txt
  # sed -n 2834p ../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/yeast_SET2del_60min_Rep2_senseStrand_full_genes_filledNAs_named.coord > 2834/2834_sense_MUT2.txt
  # sed -n 2834p ../yeast_60min_senseStrand_full_genes_filledNAs_named.coord/yeast_SET2del_60min_Rep3_senseStrand_full_genes_filledNAs_named.coord > 2834/2834_sense_MUT3.txt
  # sed -n 2834p ../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/yeast_WT_60min_Rep1_antisenseStrand_full_genes_filledNAs_named.coord > 2834/2834_antisense_WT1.txt
  # sed -n 2834p ../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/yeast_WT_60min_Rep2_antisenseStrand_full_genes_filledNAs_named.coord > 2834/2834_antisense_WT2.txt
  # sed -n 2834p ../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/yeast_WT_60min_Rep3_antisenseStrand_full_genes_filledNAs_named.coord > 2834/2834_antisense_WT3.txt
  # sed -n 2834p ../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/yeast_SET2del_60min_Rep1_antisenseStrand_full_genes_filledNAs_named.coord > 2834/2834_antisense_MUT1.txt
  # sed -n 2834p ../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/yeast_SET2del_60min_Rep2_antisenseStrand_full_genes_filledNAs_named.coord > 2834/2834_antisense_MUT2.txt
  # sed -n 2834p ../yeast_60min_antisenseStrand_full_genes_filledNAs_named.coord/yeast_SET2del_60min_Rep3_antisenseStrand_full_genes_filledNAs_named.coord > 2834/2834_antisense_MUT3.txt
  # cd 2834
  # cat *.txt > 2834.txt

  # 

### Figure 1A

singlePlot = function(Genename, linenumber, time='120min'){

x = read.table(paste0('~/Dropbox/JieHuang/plot_raw_files/', linenumber, '.txt'), header=F, sep='\t')
  
len = sum(!is.na(x[1,7:ncol(x)]))



if (time == '30min') {
  print(time)
  print('Use only WT rep1 & WT rep2')
  sense.mut = x[6:8,7:(6+len)]
  sense.wt = x[9:10,7:(6+len)]
  as.mut = x[1:3,7:(6+len)]
  as.wt = x[4:5,7:(6+len)]
} else {
  sense.mut = x[7:9,7:(6+len)]
  sense.wt = x[10:12,7:(6+len)]
  
  as.mut = x[1:3,7:(6+len)]
  as.wt = x[4:6,7:(6+len)]
}



sense.mut[1,] = ksmooth(1:len, sense.mut[1,],bandwidth = 500)$y
sense.mut[2,] = ksmooth(1:len, sense.mut[2,],bandwidth = 500)$y
sense.mut[3,] = ksmooth(1:len, sense.mut[3,],bandwidth = 500)$y
sense.wt[1,] = ksmooth(1:len, sense.wt[1,],bandwidth = 500)$y
sense.wt[2,] = ksmooth(1:len, sense.wt[2,],bandwidth = 500)$y
if (time != '30min') 
  {sense.wt[3,] = ksmooth(1:len, sense.wt[3,],bandwidth = 500)$y}



as.mut[1,] = ksmooth(1:len, as.mut[1,],bandwidth = 500)$y
as.mut[2,] = ksmooth(1:len, as.mut[2,],bandwidth = 500)$y
as.mut[3,] = ksmooth(1:len, as.mut[3,],bandwidth = 500)$y
as.wt[1,] = ksmooth(1:len, as.wt[1,],bandwidth = 500)$y
as.wt[2,] = ksmooth(1:len, as.wt[2,],bandwidth = 500)$y
if (time != '30min') 
{as.wt[3,] = ksmooth(1:len, as.wt[3,],bandwidth = 500)$y}



sense.wt.avg = apply(sense.wt, 2, mean)
sense.wt.sd = apply(sense.wt, 2, sd)
sense.mut.avg = apply(sense.mut, 2, mean)
sense.mut.sd = apply(sense.mut, 2, sd)

sense.mut.avg.ks = sense.mut.avg
sense.wt.avg.ks = sense.wt.avg
sense.mut.sd.ks = sense.mut.sd
sense.wt.sd.ks = sense.wt.sd



as.wt.avg = apply(as.wt, 2, mean)
as.wt.sd = apply(as.wt, 2, sd)
as.mut.avg = apply(as.mut, 2, mean)
as.mut.sd = apply(as.mut, 2, sd)

as.mut.avg.ks = as.mut.avg
as.wt.avg.ks = as.wt.avg
as.mut.sd.ks = as.mut.sd
as.wt.sd.ks = as.wt.sd




loss1 = read.table(paste0('~/Dropbox/JieHuang/plot_raw_files/gene_', linenumber-1, '_WTavg_MUTavg_loss.txt'), sep='\t', header =F)

loss1 = loss1 - min(loss1, na.rm=T)

#png(paste0(linenumber,'.sense.png'), width = 5*600,  height = 4*600,  res = 600, pointsize = 8)      
pdf(paste0(linenumber,'_sense_', time, '.pdf'), width = 8,  height = 12)      

#top plot for sense line plot
#par(fig=c(0.1,0.9,0.6 ,1))
par(mfrow=c(3,1))
plot(1:len, sense.mut.avg.ks,type='l',col='red', lwd=3, cex.lab=1, cex.axis=1, cex.main=1,
     ylim=c(0,80), xlim=c(0,len), 
     xlab=NA, ylab='Kernel Smoothed Read Depth (Sense)',
     main=paste0(Genename,' Kernal Smoothed Signal \n(Bandwidth 500)'))
axis(side = 2)
points(1:len, sense.wt.avg.ks , type='l',lwd=3)
polygon(c(1:len, len:1), 
        c((sense.mut.avg.ks + sense.mut.sd.ks) , rev(sense.mut.avg.ks - sense.mut.sd.ks) ),
        col=rgb(1,0,0, 0.1), border=T )
polygon(c(1:len, len:1), 
        c((sense.wt.avg.ks + sense.wt.sd.ks) , rev(sense.wt.avg.ks - sense.wt.sd.ks) ),
        col=rgb(0,0,0, 0.2), border=T )
if (linenumber==4148){abline(v=806, lty=2,col='blue')}
if (linenumber==5781){abline(v=1240, lty=2,col='blue')}
if (linenumber==2834 && time =='30min'){abline(v=1283, lty=2,col='blue')}

legend('topleft', c('WT', "SET2 deleted"), col=c('black', 'red'), lty=1, border=F, cex=0.8)

#middle plot for antisense
#par(fig=c(0.1,0.9,0.2 ,0.6))
plot(1:len, as.mut.avg.ks,type='l',col='red', lwd=3, cex.lab=1, cex.axis=1, cex.main=1,
     ylim=c(0,40), xlim=c(0,len), 
     xlab=NA, ylab='Kernel Smoothed Read Depth (Anti-Sense)')
axis(side = 2)
points(1:len, as.wt.avg.ks , type='l',lwd=3)
polygon(c(1:len, len:1), 
        c((as.mut.avg.ks + as.mut.sd.ks) , rev(as.mut.avg.ks - as.mut.sd.ks) ),
        col=rgb(1,0,0, 0.1), border=T )
polygon(c(1:len, len:1), 
        c((as.wt.avg.ks + as.wt.sd.ks) , rev(as.wt.avg.ks - as.wt.sd.ks) ),
        col=rgb(0,0,0, 0.2), border=T )

legend('topleft', c('WT', "SET2 deleted"), col=c('black', 'red'), lty=1, border=F, cex=0.8)



#bottom plot for loss band
#par(fig=c(0.1,0.9,0,0.2), new=TRUE)
plot(1:len, sense.mut.avg.ks,type='n',col='seagreen4', lwd=5, cex.lab=1, cex.axis=1, cex.main=1,
     ylim=c(0,0.1), xlim=c(0,len), yaxt='n', xaxt='n', bty="n",
     xlab='Gene (bp)', ylab='MSEL')

my_palette <- c(colorRampPalette(c("yellow", "skyblue1"))(n = 40), colorRampPalette(c("skyblue1", "blue"))(n = 70))
for (i in 151:(len-150)) {
  abline(v=i, col=my_palette[round(as.numeric(loss1[1, (i-150)]),1) * 10 +1])
}

if (linenumber==4148){abline(v=806, lty=2,col='blue')}
if (linenumber==5781){abline(v=1240, lty=2,col='blue')}

dev.off()

}

# 
# plot(1:len, c(rep(NA,149), loss1[1:(length(loss1)-1)], rep(NA, 149), loss1[length(loss1)]), 
#      ylab='Optimal MSEL', xlab='Gene (bp)',cex=1,
#      type='l', ylim=c(0, 25), xlim=c(0,len), col='forestgreen')
# abline(v=806, lty=2,col='blue')
# 
# dev.off()


singlePlot('lcb5', 4148)
singlePlot('avo1', 5781)
# singlePlot('bub1', 2695)
# singlePlot('npr2', 2970)

singlePlot('LAM5', 2834, '30min')




#png('color key.png', width = 5*600,  height = 2*600,  res = 600, pointsize = 8)    
pdf('color key.pdf', width = 5,  height = 2)    


my_palette <- c(colorRampPalette(c("yellow", "skyblue1"))(n = 40), colorRampPalette(c("skyblue1", "blue"))(n = 70))
plot((1:length(my_palette))*0.1, rep(1, length(my_palette)), 
     type='n', yaxt='n',
     xlab='MSRL', ylab='',bty="n")
for (i in 1: length(my_palette)) {
  abline(v=i*0.1, col=my_palette[i], lwd=5)
}
abline(v=7.9, col='red')
abline(v=10.4, col='red')
dev.off()





