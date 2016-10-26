#bsub Rscript ~/determine_cryptic_jumps_and_rates.R ~/single_gene_seq.txt ~/output.txt ~/loss.txt 0

# args[1] == input file, the sample file is in the same folder using gene Lcb5,
# args[2] == output file: genename, gene.len, CIS(i.e., theta), y, z, MSRLD, 
# arg[3] == output file, MSRL value of all positions except the beginning and ending 150bp, 
# args[4] == 0/1, if we want stepsize=1, args[4] = 0, if want stepsize=100, default args[4] == 1

library(lpSolve)
library(Matrix)
library(methods)


linProg = function(d0){
  len = length(d0)
  read_length = 50
  
  band.mat = matrix(0, len, len+read_length-1)
  for(j in 1:nrow(band.mat)){
    band.mat[j,(j:(j+read_length-1))] = 1
  }
  
  objective.in = rep.int(1, len+read_length-1)
  const.mat = rbind(band.mat, diag(rep.int(1, len+read_length-1)))
  const.dir = rep.int(">=", 2*len+read_length-1)
  const.rhs = c(d0, rep.int(0, len+read_length-1))
  
  return(list("band" = band.mat, "lp" = lp("min",objective.in,const.mat,const.dir,const.rhs)))
}

args = commandArgs(trailingOnly = TRUE)

avgs = read.table(args[1], sep="\t")
genename = as.character(avgs[1,1])
WT.avg = as.numeric(avgs[1,7:(ncol(avgs))])
mut.avg = as.numeric(avgs[2,7:(ncol(avgs))])

#out.mat = c("Gene_Name","gene.len", "Best_Theta","Best_y","Best_z","Max_minus_Min_of_lossFunction")
read_length = 50
  
gene.len = length(WT.avg)
print("gene length:")
print(gene.len)


if (args[4]=='1') {
  thetas = seq(151,gene.len-150,by=100)
} else if (args[4]=='0') {thetas = 151:(gene.len-150)}
thetas = c(thetas,gene.len)

par = matrix(0,length(thetas),2)
lf = matrix(0,length(thetas),1)

mut.middle = mut.avg
mut.middle[1:100] = 0
mut.middle[(gene.len-99):gene.len] = 0

linProglist = linProg(WT.avg)
band.mat = linProglist$band
dim(band.mat)
l = linProglist$lp


Xy = Matrix(band.mat, sparse = TRUE)
Xz = Matrix(band.mat, sparse = TRUE)

Xy[1:100,] = 0
Xy[(gene.len-99):gene.len,] = 0
Xz[1:100,] = 0
Xz[(gene.len-99):gene.len,] = 0

a11 = as.numeric(2*t(l$solution) %*% t(Xy) %*% Xy %*% l$solution)
b1 = as.numeric(2*t(l$solution) %*% t(Xy) %*% mut.middle)
mutmut = as.numeric(t(mut.middle) %*% mut.middle)


out.mat = c()
continue = 1
for(p in 1:(length(thetas)-1)){
  #print(thetas[p])
  Xz[,1:(thetas[p]+read_length-2)] = 0
  a12 = as.numeric(t(l$solution) %*% (t(Xy) %*% Xz + t(Xz) %*% Xy) %*% l$solution)
  a22 = as.numeric(2*t(l$solution) %*% t(Xz) %*% Xz %*% l$solution)
  b2 = as.numeric(2*t(l$solution) %*% t(Xz) %*% mut.middle)
  
  A = rbind(c(a11,a12), c(a12,a22))
  rhs = c(b1,b2)
  
  solve.ans = try(solve(A,rhs), silent = TRUE)
  if(class(solve.ans) == "try-error"){
    continue=0
    print("try error!")
    print("the first 300 lambda0:")
    print(l$solution[1:300])
    break
  }
  
  par[p,] = solve.ans
  
  if(par[p,2] < 0){
    par[p,2] = 0
    par[p,1] = b1/a11
  }
  
  y = par[p,1]
  z = par[p,2]
  #X = y*Xy + z*Xz
  lf[p] = sqrt((.5*a11*y*y + a12*y*z -b1*y + .5*a22*z*z-b2*z+mutmut)/(gene.len-200))
}

if(continue==1){
  par[length(thetas),1] = b1/a11
  #X = par[length(thetas),1]*Xy
  lf[length(thetas)] = sqrt((.5*a11*y*y + a12*y*z -b1*y + .5*a22*z*z-b2*z+mutmut)/(gene.len-200))
  
  best.theta.ind = which.min(lf)
  best.theta = thetas[best.theta.ind]
  best.y = par[best.theta.ind,1]
  best.z = par[best.theta.ind,2]

  out.mat = c(genename, gene.len, best.theta, best.y, best.z, max(lf)-min(lf))  
}else{
  out.mat = c(genename, "Singular", "NA", "NA", "NA", "NA")
}

write.table(as.matrix(t(out.mat)), file=args[2], sep="\t", col.names=F, row.names=F, quote=F)

write.table(as.matrix(t(lf)), file=args[3], sep="\t", col.names=F, row.names=F, quote=F)
