mutPlot = function(time){
  print(time)
  ### sense, mut
  if (time == '0min'){
    mut2 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep2_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut3 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep3_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut4_inf = rbind( mut2[1,], mut3[1,])
    mut2_4 = rbind( mut2[2,], mut3[2,])
    mut0_2 = rbind( mut2[3,], mut3[3,])
  } else if (time=='0min_WT1_mut2'){
    mut2 = read.table('~/Downloads/yeast_SET2del_0min_Rep2_senseStrand_full_genes_filledNAs_named.coord_average_WT1_mut2.txt', header=T, sep='\t')
    mut4_inf = mut2[1,]
    mut2_4 = mut2[2,]
    mut0_2 = mut2[3,]
  }else {
    mut1 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep1_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut2 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep2_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut3 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep3_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut4_inf = rbind(mut1[1,], mut2[1,], mut3[1,])
    mut2_4 = rbind(mut1[2,], mut2[2,], mut3[2,])
    mut0_2 = rbind(mut1[3,], mut2[3,], mut3[3,])
  }
  
  mut4_inf_avg = apply(mut4_inf,2,mean)
  mut2_4_avg = apply(mut2_4,2,mean)
  mut0_2_avg = apply(mut0_2,2,mean)
  
  mut4_inf_sd = apply(mut4_inf, 2, sd)
  mut2_4_sd = apply(mut2_4,2,sd)
  mut0_2_sd = apply(mut0_2,2,sd)
  
  mut4_inf_avg = mut4_inf_avg - min(mut4_inf_avg, na.rm=T)
  mut2_4_avg = mut2_4_avg - min(mut2_4_avg, na.rm=T)
  mut0_2_avg = mut0_2_avg - min(mut0_2_avg, na.rm=T)
  
  
  
  #antisense, mut
  
  if (time == '0min'){
    mut2.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep2_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut3.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep3_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut4_inf.as = rbind( mut2.as[1,], mut3.as[1,])
    mut2_4.as = rbind( mut2.as[2,], mut3.as[2,])
    mut0_2.as = rbind( mut2.as[3,], mut3.as[3,])
    
  } else if (time=='0min_WT1_mut2'){
    mut2.as = read.table('~/Downloads/yeast_SET2del_0min_Rep2_antisenseStrand_full_genes_filledNAs_named.coord_average_0min_WT1_mut2.txt', header=T, sep='\t')
    mut4_inf.as = mut2.as[1,]
    mut2_4.as = mut2.as[2,]
    mut0_2.as = mut2.as[3,]
  } else {
    mut1.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep1_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut2.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep2_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    mut3.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_SET2del_', time, '_Rep3_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    mut4_inf.as = rbind(mut1.as[1,], mut2.as[1,], mut3.as[1,])
    mut2_4.as = rbind(mut1.as[2,], mut2.as[2,], mut3.as[2,])
    mut0_2.as = rbind(mut1.as[3,], mut2.as[3,], mut3.as[3,])
  }
  
  
  
  
  mut4_inf_avg.as = apply(mut4_inf.as,2,mean)
  mut2_4_avg.as = apply(mut2_4.as,2,mean)
  mut0_2_avg.as = apply(mut0_2.as,2,mean)
  
  mut4_inf_sd.as = apply(mut4_inf.as, 2, sd)
  mut2_4_sd.as = apply(mut2_4.as,2,sd)
  mut0_2_sd.as = apply(mut0_2.as,2,sd)
  
  mut4_inf_avg.as = mut4_inf_avg.as - min(mut4_inf_avg.as, na.rm=T)
  mut2_4_avg.as = mut2_4_avg.as - min(mut2_4_avg.as, na.rm=T)
  mut0_2_avg.as = mut0_2_avg.as - min(mut0_2_avg.as, na.rm=T)
  
  
  J=1000; K=1000;
  library(scales)
  #png(paste0('~/Downloads/', time, '_mut.png'), width = 5*600,  height = 5*600,  res = 600, pointsize = 8)   
  pdf(paste0('~/Downloads/', time, '_mut.pdf'), width = 6,  height = 8)   
  par(mfrow=c(2,1))
  plot(-J:K, mut0_2_avg, type='l',col='red',ylim=c(0,40), xlim=c(-J,K), lwd=2,
       ylab='Baseline Normalized Read Depth', xlab='Position to Jump Point')
  abline(v=0,lty=2)
  polygon(c(-J:K, rev(-J:K)), c(mut0_2_avg+mut0_2_sd, rev(mut0_2_avg-mut0_2_sd)),col = alpha('red', 0.3), border = F)
  
  polygon(c((-J:K)[!is.na(mut2_4_avg)], rev((-J:K)[!is.na(mut2_4_avg)])), c((mut2_4_avg+mut2_4_sd)[!is.na(mut2_4_avg)], rev((mut2_4_avg-mut2_4_sd)[!is.na(mut2_4_avg)])),col = alpha('black', 0.3), border = F)
  lines(-J:K, mut2_4_avg, type='l',col='black', lwd=2)
  
  polygon(c((-J:K)[!is.na(mut4_inf_avg)], rev((-J:K)[!is.na(mut4_inf_avg)])), c((mut4_inf_avg+mut4_inf_sd)[!is.na(mut4_inf_avg)], rev((mut4_inf_avg-mut4_inf_sd)[!is.na(mut4_inf_avg)])),col = alpha('blue', 0.3), border = F)
  lines(-J:K, mut4_inf_avg, type='l',col='blue', lwd=2)
  
  legend('topleft', c('High','Medium','Low'), col = c('blue','black','red'),lty = c(1, 1, 1), title='Cryp Init Chance', cex=0.7)
  
  #text(400, 50, "Sense Strand", cex=1.5)
  
  title(paste0('SET2 deleted group\n ', time, ' after nutrition deprivation'))
  
  
  #_______antisense plot
  plot(-J:K, mut0_2_avg.as, type='l',col=alpha('red', 0.8), ylim=c(0,15), xlim=c(-J,K), lwd=2,
       ylab='Baseline Normalized Read Depth', xlab='Position to Jump Point')
  polygon(c(-J:K, rev(-J:K)), c(mut0_2_avg.as+mut0_2_sd.as, rev(mut0_2_avg.as-mut0_2_sd.as)),col = alpha('red', 0.1), border = F)
  
  
  polygon(c((-J:K)[!is.na(mut2_4_avg.as)], rev((-J:K)[!is.na(mut2_4_avg.as)])), c((mut2_4_avg.as+mut2_4_sd.as)[!is.na(mut2_4_avg.as)], rev((mut2_4_avg.as-mut2_4_sd.as)[!is.na(mut2_4_avg.as)])),col = alpha('black', 0.1), border = F)
  lines(-J:K, mut2_4_avg.as, type='l',col=alpha('black',0.8), lwd=2)
  
  polygon(c((-J:K)[!is.na(mut4_inf_avg.as)], rev((-J:K)[!is.na(mut4_inf_avg.as)])), c((mut4_inf_avg.as+mut4_inf_sd.as)[!is.na(mut4_inf_avg.as)], rev((mut4_inf_avg.as-mut4_inf_sd.as)[!is.na(mut4_inf_avg.as)])),col = alpha('blue', 0.1), border = F)
  lines(-J:K, mut4_inf_avg.as, type='l',col=alpha('blue',0.8), lwd=2)
  
  abline(v=0,lty=2)
  
  #text(600, 15, "Antisense Strand", cex=1.5)
  
  
  #legend('topleft', c('High','Medium','Low'), col = c('blue','black','red'),lty = c(1, 1, 1), title='Cryp Init Chance', cex=0.7)
  
  dev.off()
  
  
}

##################################################

WTPlot = function(time){
  print(time)
  if (time == '0min'){
    WT3 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep3_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT4 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep4_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT5 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep5_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT6 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep6_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    WT4_inf = rbind( WT3[1,], WT4[1,], WT5[1,], WT6[1,])
    WT2_4 = rbind( WT3[2,], WT4[2,], WT5[2,], WT6[2,])
    WT0_2 = rbind( WT3[3,], WT4[3,], WT5[3,], WT6[3,])
  } else if (time=='0min_WT1_mut2'){
    WT2 = read.table('~/Downloads/yeast_WT_0min_Rep1_senseStrand_full_genes_filledNAs_named.coord_average_WT1_mut2.txt', header=T, sep='\t')
    WT4_inf = WT2[1,]
    WT2_4 = WT2[2,]
    WT0_2 = WT2[3,]
  } else if (time == '30min'){
    WT1 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep1_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT2 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep2_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    WT4_inf = rbind(WT1[1,], WT2[1,])
    WT2_4 = rbind(WT1[2,], WT2[2,])
    WT0_2 = rbind(WT1[3,], WT2[3,])    
  } else {
    WT1 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep1_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT2 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep2_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT3 = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep3_senseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    WT4_inf = rbind(WT1[1,], WT2[1,], WT3[1,])
    WT2_4 = rbind(WT1[2,], WT2[2,], WT3[2,])
    WT0_2 = rbind(WT1[3,], WT2[3,], WT3[3,])
  }
  
  WT4_inf_avg = apply(WT4_inf,2,mean)
  WT2_4_avg = apply(WT2_4,2,mean)
  WT0_2_avg = apply(WT0_2,2,mean)
  
  WT4_inf_sd = apply(WT4_inf, 2, sd)
  WT2_4_sd = apply(WT2_4,2,sd)
  WT0_2_sd = apply(WT0_2,2,sd)
  
  WT4_inf_avg = WT4_inf_avg - min(WT4_inf_avg, na.rm=T)
  WT2_4_avg = WT2_4_avg - min(WT2_4_avg, na.rm=T)
  WT0_2_avg = WT0_2_avg - min(WT0_2_avg, na.rm=T)
  
  
  
  #antisense, WT
  
  if (time == '0min'){
    WT3.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep3_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT4.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep4_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT5.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep5_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT6.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep6_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    WT4_inf.as = rbind( WT3.as[1,], WT4.as[1,],WT5.as[1,], WT6.as[1,])
    WT2_4.as = rbind( WT3.as[2,], WT4.as[2,],WT5.as[2,], WT6.as[2,])
    WT0_2.as = rbind( WT3.as[3,], WT4.as[3,],WT5.as[3,], WT6.as[3,])
    
  } else if (time=='0min_WT1_mut2'){
    WT2.as = read.table('~/Downloads/yeast_WT_0min_Rep1_antisenseStrand_full_genes_filledNAs_named.coord_average_WT1_mut2.txt', header=T, sep='\t')
    WT4_inf.as = WT2.as[1,]
    WT2_4.as = WT2.as[2,]
    WT0_2.as = WT2.as[3,]
  } else if (time == '30min'){
    WT1.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep1_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT2.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep2_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    WT4_inf.as = rbind(WT1.as[1,], WT2.as[1,])
    WT2_4.as = rbind(WT1.as[2,], WT2.as[2,])
    WT0_2.as = rbind(WT1.as[3,], WT2.as[3,])
  } else {
    WT1.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep1_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT2.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep2_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    WT3.as = read.table(paste0('/Users/JieHuang1/Dropbox/JieHuang/plot_raw_files/Panel_D_lineplot/fix_th_yeast_WT_', time, '_Rep3_antisenseStrand_full_genes_filledNAs_named.coord_upstream1000_downstream1000_average.txt'), header=T, sep='\t')
    
    WT4_inf.as = rbind(WT1.as[1,], WT2.as[1,], WT3.as[1,])
    WT2_4.as = rbind(WT1.as[2,], WT2.as[2,], WT3.as[2,])
    WT0_2.as = rbind(WT1.as[3,], WT2.as[3,], WT3.as[3,])
  }
  
  
  
  
  WT4_inf_avg.as = apply(WT4_inf.as,2,mean)
  WT2_4_avg.as = apply(WT2_4.as,2,mean)
  WT0_2_avg.as = apply(WT0_2.as,2,mean)
  
  WT4_inf_sd.as = apply(WT4_inf.as, 2, sd)
  WT2_4_sd.as = apply(WT2_4.as,2,sd)
  WT0_2_sd.as = apply(WT0_2.as,2,sd)
  
  WT4_inf_avg.as = WT4_inf_avg.as - min(WT4_inf_avg.as, na.rm=T)
  WT2_4_avg.as = WT2_4_avg.as - min(WT2_4_avg.as, na.rm=T)
  WT0_2_avg.as = WT0_2_avg.as - min(WT0_2_avg.as, na.rm=T)
  
  
  J=1000; K=1000;
  library(scales)
  #png(paste0('~/Downloads/', time, '_WT.png'), width = 5*600,  height = 5*600,  res = 600, pointsize = 8)   
  pdf(paste0('~/Downloads/', time, '_WT.pdf'), width = 6,  height = 8)  
  par(mfrow=c(2,1))
  plot(-J:K, WT0_2_avg, type='l',col='red',ylim=c(0,40), xlim=c(-J,K), lwd=2,
       ylab='Baseline Normalized Read Depth', xlab='Position to Jump Point')
  abline(v=0,lty=2)
  polygon(c(-J:K, rev(-J:K)), c(WT0_2_avg+WT0_2_sd, rev(WT0_2_avg-WT0_2_sd)),col = alpha('red', 0.3), border = F)
  
  polygon(c((-J:K)[!is.na(WT2_4_avg)], rev((-J:K)[!is.na(WT2_4_avg)])), c((WT2_4_avg+WT2_4_sd)[!is.na(WT2_4_avg)], rev((WT2_4_avg-WT2_4_sd)[!is.na(WT2_4_avg)])),col = alpha('black', 0.3), border = F)
  lines(-J:K, WT2_4_avg, type='l',col='black', lwd=2)
  
  polygon(c((-J:K)[!is.na(WT4_inf_avg)], rev((-J:K)[!is.na(WT4_inf_avg)])), c((WT4_inf_avg+WT4_inf_sd)[!is.na(WT4_inf_avg)], rev((WT4_inf_avg-WT4_inf_sd)[!is.na(WT4_inf_avg)])),col = alpha('blue', 0.3), border = F)
  lines(-J:K, WT4_inf_avg, type='l',col='blue', lwd=2)
  
  legend('topleft', c('High','Medium','Low'), col = c('blue','black','red'),lty = c(1, 1, 1), title='Cryp Init Chance', cex=0.7)
  
  #text(600, 50, "Sense Strand", cex=1.5)
  
  title(paste0('Wild Type group\n ', time, ' after nutrition deprivation'))
  
  
  #_______antisense plot
  plot(-J:K, WT0_2_avg.as, type='l',col=alpha('red', 0.8), ylim=c(0,15), xlim=c(-J,K), lwd=2,
       ylab='Baseline Normalized Read Depth', xlab='Position to Jump Point')
  polygon(c(-J:K, rev(-J:K)), c(WT0_2_avg.as+WT0_2_sd.as, rev(WT0_2_avg.as-WT0_2_sd.as)),col = alpha('red', 0.1), border = F)
  
  
  polygon(c((-J:K)[!is.na(WT2_4_avg.as)], rev((-J:K)[!is.na(WT2_4_avg.as)])), c((WT2_4_avg.as+WT2_4_sd.as)[!is.na(WT2_4_avg.as)], rev((WT2_4_avg.as-WT2_4_sd.as)[!is.na(WT2_4_avg.as)])),col = alpha('black', 0.1), border = F)
  lines(-J:K, WT2_4_avg.as, type='l',col=alpha('black',0.8), lwd=2)
  
  polygon(c((-J:K)[!is.na(WT4_inf_avg.as)], rev((-J:K)[!is.na(WT4_inf_avg.as)])), c((WT4_inf_avg.as+WT4_inf_sd.as)[!is.na(WT4_inf_avg.as)], rev((WT4_inf_avg.as-WT4_inf_sd.as)[!is.na(WT4_inf_avg.as)])),col = alpha('blue', 0.1), border = F)
  lines(-J:K, WT4_inf_avg.as, type='l',col=alpha('blue',0.8), lwd=2)
  
  abline(v=0,lty=2)
  
  #text(600, 15, "Antisense Strand", cex=1.5)
  
  
  #legend('topleft', c('High','Medium','Low'), col = c('blue','black','red'),lty = c(1, 1, 1), title='Cryp Init Chance', cex=0.7)
  
  dev.off()
  
  
}



mutPlot('0min')
mutPlot('30min')
mutPlot('60min')
mutPlot('120min')
# mutPlot('0min_WT1_mut2')


WTPlot('0min')
WTPlot('30min')
WTPlot('60min')
WTPlot('120min')
# WTPlot('0min_WT1_mut2')





