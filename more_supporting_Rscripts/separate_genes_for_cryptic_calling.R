# for 0 mins, using WT rep 3,4,5,6, SET2del rep2,3

# args[1] == WT_rep3, args[2] == WT_rep4, args[3] == WT_rep5,args[4] == WT_rep6, 
# args[5] == Mut_rep2, args[6] == Mut_rep3, 
# args[7] == "0min"

# submit code
# 0 min
#cd /proj/strahllb/users/Jie/Stephen/yeast_0min_senseStrand_full_genes_filledNAs_named.coord
#bsub -M 8 -q week Rscript /proj/strahllb/users/Jie/Stephen/separate_genes_for_cryptic_calling.R  yeast_WT_0min_Rep3_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_0min_Rep4_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_0min_Rep5_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_0min_Rep6_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_0min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_0min_Rep3_senseStrand_full_genes_filledNAs_named.coord 0min


# for 30min
# args[1] == WT_rep1, args[2] == WT_rep2, args[3] == WT_rep3 (not using WT_rep3),
# args[3] == Mut_rep1,args[4] == Mut_rep2,args[5] == Mut_rep3, 
# args[7] == "30min"
#cd /proj/strahllb/users/Jie/Stephen/yeast_30min_senseStrand_full_genes_filledNAs_named.coord 
#bsub -M 8 -q week Rscript /proj/strahllb/users/Jie/Stephen/separate_genes_for_cryptic_calling.R  yeast_WT_30min_Rep1_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_30min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_30min_Rep3_senseStrand_full_genes_filledNAs_named.coord yeast_SET2del_30min_Rep1_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_30min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_30min_Rep3_senseStrand_full_genes_filledNAs_named.coord 30min 

 

# for  60min, 120min, WT rep1,2,3, SETdel rep1,2,3
# args[1] == WT_rep1, args[2] == WT_rep2, args[3] == WT_rep3,
# args[4] == Mut_rep1,args[5] == Mut_rep2,args[6] == Mut_rep3, 
# args[7] == "60min" 

# submit code
# for 60min
#cd /proj/strahllb/users/Jie/Stephen/yeast_60min_senseStrand_full_genes_filledNAs_named.coord 
#bsub -M 8 -q week Rscript /proj/strahllb/users/Jie/Stephen/separate_genes_for_cryptic_calling.R  yeast_WT_60min_Rep1_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_60min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_60min_Rep3_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_60min_Rep1_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_60min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_60min_Rep3_senseStrand_full_genes_filledNAs_named.coord 60min 

# for 120min
#cd /proj/strahllb/users/Jie/Stephen/yeast_120min_senseStrand_full_genes_filledNAs_named.coord 
#bsub -M 8 -q week Rscript /proj/strahllb/users/Jie/Stephen/separate_genes_for_cryptic_calling.R  yeast_WT_120min_Rep1_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_120min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_WT_120min_Rep3_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_120min_Rep1_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_120min_Rep2_senseStrand_full_genes_filledNAs_named.coord  yeast_SET2del_120min_Rep3_senseStrand_full_genes_filledNAs_named.coord 120min 



# Jan 23 2017
# Using data from paper: Selective suppression of antisense transcription by Set2-mediated H3K36 methylation, Workman, 2016
# args[1] == WT_rep1, args[2] == WT_rep2, args[3] = not using
# args[4] == Mut_rep1,args[5] == Mut_rep2,args[6] == not using, 
# args[7] == "Workman"
# cd /proj/strahllb/users/Jie/Stephen/Workman2016/
# bsub -M 8 -q week Rscript /proj/strahllb/users/Jie/Stephen/Rscripts/separate_genes_for_cryptic_calling.R  \
# /proj/strahllb/users/Austin/Stephen/nutrient_deprivation/seq_removedDups/SRAT_paper/bowtie_out/sense/coords/yeast_BY4741_WT_Rep1_full_gene_senseStrand_filledNAs.coord \
# /proj/strahllb/users/Austin/Stephen/nutrient_deprivation/seq_removedDups/SRAT_paper/bowtie_out/sense/coords/yeast_BY4741_WT_Rep2_full_gene_senseStrand_filledNAs.coord \
# /proj/strahllb/users/Austin/Stephen/nutrient_deprivation/seq_removedDups/SRAT_paper/bowtie_out/sense/coords/yeast_BY4741_WT_Rep2_full_gene_senseStrand_filledNAs.coord \
# /proj/strahllb/users/Austin/Stephen/nutrient_deprivation/seq_removedDups/SRAT_paper/bowtie_out/sense/coords/yeast_BY4741_SET2del_Rep1_full_gene_senseStrand_filledNAs.coord \
# /proj/strahllb/users/Austin/Stephen/nutrient_deprivation/seq_removedDups/SRAT_paper/bowtie_out/sense/coords/yeast_BY4741_SET2del_Rep2_full_gene_senseStrand_filledNAs.coord \
# /proj/strahllb/users/Austin/Stephen/nutrient_deprivation/seq_removedDups/SRAT_paper/bowtie_out/sense/coords/yeast_BY4741_SET2del_Rep2_full_gene_senseStrand_filledNAs.coord \
# Workman   

#set.tempdir('/netscr/jh259/Workman2016') 

args = commandArgs(trailingOnly = TRUE)

#read in matrices
if (args[7] == "0min") {
    print("Data from:")
    print(args[7])
	WT_rep1.mat = read.table(args[1], sep="\t", row.names=1, header=T)
	WT_rep2.mat = read.table(args[2], sep="\t", row.names=1, header=T)
	WT_rep3.mat = read.table(args[3], sep="\t", row.names=1, header=T)
	WT_rep4.mat = read.table(args[4], sep="\t", row.names=1, header=T)
	mut_rep1.mat = read.table(args[5], sep="\t", row.names=1, header=T)
	mut_rep2.mat = read.table(args[6], sep="\t", row.names=1, header=T)
} else if (args[7] == "30min") {
    print("Data from:")
    print(args[7])
    WT_rep1.mat = read.table(args[1], sep="\t", row.names=1, header=T)
	WT_rep2.mat = read.table(args[2], sep="\t", row.names=1, header=T)
	mut_rep1.mat = read.table(args[4], sep="\t", row.names=1, header=T)
	mut_rep2.mat = read.table(args[5], sep="\t", row.names=1, header=T)
	mut_rep3.mat = read.table(args[6], sep="\t", row.names=1, header=T)
} else if (args[7] == "Workman") {
    print("Data from:")
    print(args[7])
    WT_rep1.mat = read.table(args[1], sep="\t", header=T)
	WT_rep2.mat = read.table(args[2], sep="\t", header=T)
	mut_rep1.mat = read.table(args[4], sep="\t", header=T)
	mut_rep2.mat = read.table(args[5], sep="\t", header=T)
	ttt = read.table('/proj/strahllb/users/Jie/Stephen/Workman2016/Gene_names.txt', sep='\t', header=F)
	row.names(WT_rep1.mat)=as.character(as.list(t(ttt)))
	row.names(WT_rep2.mat)=as.character(as.list(t(ttt)))
	row.names(mut_rep1.mat)=as.character(as.list(t(ttt)))
	row.names(mut_rep2.mat)=as.character(as.list(t(ttt)))
} else {
    print("Data from:")
    print(args[7])
    WT_rep1.mat = read.table(args[1], sep="\t", row.names=1, header=T)
	WT_rep2.mat = read.table(args[2], sep="\t", row.names=1, header=T)
	WT_rep3.mat = read.table(args[3], sep="\t", row.names=1, header=T)
	mut_rep1.mat = read.table(args[4], sep="\t", row.names=1, header=T)
	mut_rep2.mat = read.table(args[5], sep="\t", row.names=1, header=T)
	mut_rep3.mat = read.table(args[6], sep="\t", row.names=1, header=T)
}

print(WT_rep1.mat[1:5,1:10])

overlap = read.table('/proj/strahllb/users/Jie/Stephen/overlapped_gene_names.txt', sep='\t')
intron = read.table('/proj/strahllb/users/Jie/Stephen/intron_gene_names.txt', sep='\t')
chrM = read.table('/proj/strahllb/users/Jie/Stephen/chrM_gene_names.txt', sep='\t')


out.mat = c()

#for(i in 1:nrow(WT_rep1.mat)){
for(i in 5000:6000){
   #ignore overlap, intron and chrM genes:

  if (row.names(WT_rep1.mat)[i] %in% overlap[,2]) {
   cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " is a overlapped gene. \n"))
   out.mat = rbind(out.mat, c(row.names(WT_rep1.mat)[i],  "overlap ","NA", "NA", "NA", "NA"))
   next
   }
  else if (row.names(WT_rep1.mat)[i] %in% intron[,2]){
   cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " has intron(s) in it. \n"))
   out.mat = rbind(out.mat, c(row.names(WT_rep1.mat)[i],  "has intron ","NA", "NA", "NA", "NA"))
   next
   } 
  else if (row.names(WT_rep1.mat)[i] %in% chrM[,2]) {
   cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " is a mitochondrion gene. \n"))
   out.mat = rbind(out.mat, c(row.names(WT_rep1.mat)[i],  "on chrM ","NA", "NA", "NA", "NA"))
   next
   }
   
  #get average of reps
	if (args[7] == "0min") {
	  WT.avg = apply(rbind(WT_rep1.mat[i,(6:ncol(WT_rep1.mat))],WT_rep2.mat[i,(6:ncol(WT_rep2.mat))],WT_rep3.mat[i,(6:ncol(WT_rep3.mat))],WT_rep4.mat[i,(6:ncol(WT_rep4.mat))]),2, mean)
	  WT.avg = WT.avg[!is.na(WT.avg)]
	  mut.avg = apply(rbind(mut_rep1.mat[i,(6:ncol(mut_rep1.mat))],mut_rep2.mat[i,(6:ncol(mut_rep2.mat))]),2, mean)
	  mut.avg = mut.avg[!is.na(mut.avg)]
	} else if (args[7] == "30min") {
	  WT.avg = apply(rbind(WT_rep1.mat[i,(6:ncol(WT_rep1.mat))],WT_rep2.mat[i,(6:ncol(WT_rep2.mat))]),2, mean)
	  WT.avg = WT.avg[!is.na(WT.avg)]
	  mut.avg = apply(rbind(mut_rep1.mat[i,(6:ncol(mut_rep1.mat))],mut_rep2.mat[i,(6:ncol(mut_rep2.mat))],mut_rep3.mat[i,(6:ncol(mut_rep3.mat))]),2, mean)
	  mut.avg = mut.avg[!is.na(mut.avg)]
	} else if (args[7] == "Workman") {
	  WT.avg = apply(rbind(WT_rep1.mat[i,(7:ncol(WT_rep1.mat))],WT_rep2.mat[i,(7:ncol(WT_rep2.mat))]),2, mean)
	  WT.avg = WT.avg[!is.na(WT.avg)]
	  mut.avg = apply(rbind(mut_rep1.mat[i,(7:ncol(mut_rep1.mat))],mut_rep2.mat[i,(7:ncol(mut_rep2.mat))]),2, mean)
	  mut.avg = mut.avg[!is.na(mut.avg)]
	} else {
	  WT.avg = apply(rbind(WT_rep1.mat[i,(6:ncol(WT_rep1.mat))],WT_rep2.mat[i,(6:ncol(WT_rep2.mat))],WT_rep3.mat[i,(6:ncol(WT_rep3.mat))]),2, mean)
	  WT.avg = WT.avg[!is.na(WT.avg)]
	  mut.avg = apply(rbind(mut_rep1.mat[i,(6:ncol(mut_rep1.mat))],mut_rep2.mat[i,(6:ncol(mut_rep2.mat))],mut_rep3.mat[i,(6:ncol(mut_rep3.mat))]),2, mean)
	  mut.avg = mut.avg[!is.na(mut.avg)]
	}
  
  if(mean(WT.avg) <= 5 | mean(mut.avg) <= 5){
    cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " had too few reads! It has an average of ", mean(WT.avg), " WT reads and ", mean(mut.avg), " mutant reads; it must be at least 5\n"))
    out.mat = rbind(out.mat, c(row.names(WT_rep1.mat)[i],  "reads<=5 ","NA", "NA", "NA", "NA"))
    next
  }
  else if(length(WT.avg) <= 700){
    cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " was too small! It was ", length(WT.avg), " and it must be at least 700bp\n"))
    out.mat = rbind(out.mat, c(row.names(WT_rep1.mat)[i], "lens<=700 ", "NA", "NA", "NA", "NA"))
    next
  }
  
  else if(length(WT.avg) > 5000){
    cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " length > 5000! It was ", length(WT.avg), "\n"))
    write.table(rbind(c(row.names(WT_rep1.mat[i,]), WT_rep1.mat[i,1:6], WT.avg), c(row.names(WT_rep1.mat[i,]), WT_rep1.mat[i,1:6], mut.avg)), file=paste0("gene_",i,"_WTavg_MUTavg.txt"), sep="\t", quote=F, row.names=F, col.names=F)
    system(paste0("bsub -M 20 -q week -J gene",i," Rscript /proj/strahllb/users/Jie/Stephen/Rscripts/determine_cryptic_jumps_and_rates.R gene_",i,"_WTavg_MUTavg.txt gene_",i,"_WTavg_MUTavg_out.txt gene_",i,"_WTavg_MUTavg_loss.txt 0"))
    next
  }
  
  cat(paste0("Gene #", i, " aka ", row.names(WT_rep1.mat)[i], " is the ones we would study. \n"))
  write.table(rbind(c(row.names(WT_rep1.mat[i,]), WT_rep1.mat[i,1:6], WT.avg), c(row.names(WT_rep1.mat[i,]), WT_rep1.mat[i,1:6], mut.avg)), file=paste0("gene_",i,"_WTavg_MUTavg.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  system(paste0("bsub -M 4 -q week -J gene",i," Rscript /proj/strahllb/users/Jie/Stephen/Rscripts/determine_cryptic_jumps_and_rates.R gene_",i,"_WTavg_MUTavg.txt gene_",i,"_WTavg_MUTavg_out.txt gene_",i,"_WTavg_MUTavg_loss.txt 0"))
}

colnames(out.mat)= c("Gene_Name","notes", "Best_Theta","Best_y","Best_z","Max_minus_Min_of_lossFunction")


#write.table(as.matrix(out.mat), file="failed_genes_WTavg_MUTavg_out.txt", quote=F, sep="\t", col.names=T, row.names=F)
write.table(as.matrix(out.mat), file="failed_genes_WTavg_MUTavg_partial_out.txt", quote=F, sep="\t", col.names=T, row.names=F)
