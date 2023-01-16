rm(list=ls())

# This repository contains the analysis pipeline of detecting marginally strong and marginally weak signals for binary outcome classification using single nucleotide polymorphism (SNP) genotype data. The developed analytical pipeline is open-source, and can be used for analyzing SNP data in R softwares with step-by-step instructions.


# Overview of Workflow and Pipeline
# 0) Housekeeping
# 1) Load outcome labels and demographic data of study population
# 2) Filter SNP genotype data using minor allele frequency threshold
# 3) Select marginally strong SNPs using the filtered genotype data
# 4) Search for marginally strong SNP correlated SNPs based on LD
# 5) Construct SNP clusters based on LD matrices
# 6) Prepare data which contains all SNPs in the detected SNP clusters 
# 7) Detect marginally weak signals based on cluster-adjusted effect sizes
# 8) Binary outcome classification based on selected SNPs
# 9) Permutation p value for each individual SNP


# Shorthands: 
# Single nucleotide polymorphism (SNP)
# Chromosome (chr)
# Minor allele frequency (MAF)
# Linkage disequilibrium (LD)

 
 


# 0 Housekeeping  --------------------------------------------------------------

## 0.1 R packages will be used -----
library(snpStats) 
library(dplyr)
library(stringr)
library(Matrix)
library(matrixcalc)
library(ggplot2)

## 0.2 User input dataset ------
### Folder path
dir.data <- file.path("/home/xuz37/AMD/34Regions/Data") 
dir.output <- file.path("/home/xuz37/AMD/34Regions/Data/Pipeline")

### Genotype Plink file name, which contains .bed, .bim, .fam  files for each chr 
bfilename = "IAMDGC" 

## 0.3 Parameters ---------
# MAF threshold 
mafthre = 0.10

# Threshold of marginal test p value for selecting marginally strong SNPs
pvalue.strong.thred = 5e-08  

# Threshold of LD, which is used to extract SNPs connected with the marginally strong SNPs
LD.thre = 0.5 

# Threshold of LD, which is used to further filter SNPs connected with marginally strong SNPs (if needed)
alpha = 0.7






# 1. Load outcome labels and demographic data of study population -------------------------
# In this section, we prepare the outcome labels and demographic data for our analysis. A vector called 'Label' need to be created to store the binary outcome data. Demographic data can be used if available. If there is any missing data in outcome labels or demographic information, please choose appropriate missing data imputation method to impute the data. Besides, we need to make the orders of subjects in labels, demographic and genotype data the same. 



## 1.1 Load demographic data (called 'Demog'), which contains  -----------------
# Label contains the binary outcome data. 
# Demographics info. For example, here we use 'age' and 'sex'. 
 
## 1.2 Make sure the orders of subjects are matched in the demographic and genotype data sets ----
# Subjects from the genotype files
## Note: Use a chr with smaller number of SNPs to save time and memory. 
setwd(dir.data)
tmp.chr = 22
tmp.genotype <- snpStats::read.plink(
                   bed=paste0("MAF", mafthre,"_CHR" ,tmp.chr,".bed"),
                   bim=paste0("MAF", mafthre,"_CHR", tmp.chr,".bim"),
                   fam=paste0("MAF", mafthre,"_CHR", tmp.chr,".fam"),
                   na.strings = c("0", "-9"), sep = "." ,
                   select.subjects = NULL, select.snps = NULL)
Genotype_ID  <- row.names(as.data.frame(as.matrix(tmp$genotypes)))
rm(list = ls(pattern="^tmp"))

# Make the orders of subjects matched in the two data sets 
Demog <- Demog[order(match(Demog$SUBJECT_ID, Genotype_ID)), ] 

# Subject ID (in order)  
subject.ID <- Demog$SUBJECT_ID

## 1.3 Outcome of interest -----
Label <- Demog$AFFECTION_STATUS 





# 2. Filter genotype data using MAF threshold----------------------------------
# In this section, we use PLINK (trough R) to filter the genotype data. The filtering is done for each chromosome separately. Basically, the SNPs with MAF larger than the pre-defined MAF threshold will be selected. The selected SNPs will be stored into a new file. 


# A list to store the number of SNPs with MAF > MAF threshold on each chr
bed.dim.list<- list() 

# Filtering data 
for(chr in 1:22){ print(chr)
  # Filtering data by calling PLINK in R
  setwd(dir.data)
  cmd <- str_c("plink --bfile ", bfilename, "_CHR", chr,"_c12_white.dose ", # Original filename 
               "--keep SubjectID.txt ",
               "--maf " , mafthre, " ",  # MAF threhod 
               "--make-bed --out MAF", mafthre, "_CHR", chr # output filename 
  )
  system(cmd)
  
  # Read in filtered data to get a sense of the dimension of filtered data by chr
  tmp <- snpStats::read.plink( 
                     bed=paste0("MAF", mafthre,"_CHR" ,chr,".bed"),
                     bim=paste0("MAF", mafthre,"_CHR", chr,".bim"),
                     fam=paste0("MAF", mafthre,"_CHR", chr,".fam"),
                     na.strings = c("0", "-9"), sep = "." , 
                     select.subjects = NULL, select.snps = NULL)
  tmp.gene <- as.data.frame(as.matrix(tmp$genotypes))
  bed.dim.list[[chr]] <-  dim(tmp.gene)
}

bed.dim.list # List, each element contains the dimension of filtered data in each chr 



# 3. Select marginally strong SNPs using the filtered genotype data-------------
# In this section, marginal test is conducted for each SNP with MAF larger than MAF threshold. Then, based on the marginally p-values and the pre-determined p-value threshold, marginally strong SNPs are determined. More specifically, genotype data are treated under additive genetic model. R package 'snpStats' is used to read the SNP information into R, which 00 represents missing data and is imputed for each SNP. Marginally test statistics and associated p-values are stored in .csv files.

## 3.1 Marginal tests --------
# There are many SNPs in the first 6 chrs. To save memory, I splitted each of these chr. into smaller working blocks 
for(chr in 1:22 ){  
  if (chr < 7){  
    tmp.size <- bed.dim.list[[chr]][2]
    l <- split(1:tmp.size, ceiling(seq_along(1:tmp.size)/100000))
    for(l.idx in 1:length(l)){
      selectsnp <- l[[l.idx]]
      # -------------------------------------------- load in data
      setwd(dir.data)
      tmp <- snpStats::read.plink(
                         bed=paste0("MAF", mafthre,"_CHR" ,chr,".bed"),
                         bim=paste0("MAF", mafthre,"_CHR", chr,".bim"),
                         fam=paste0("MAF", mafthre,"_CHR", chr,".fam"),
                         na.strings = c("0", "-9"), sep = "." ,
                         select.subjects = NULL, select.snps = selectsnp)
      tmp.gene <- as.data.frame(as.matrix(tmp$genotypes))  
      tmp.gene <- tmp.gene[order(match(rownames(tmp.gene), subject.ID)), ]
      rm(list=c("tmp", "tmp.gene"))
      
      # -------------------------------------------- Missing genotype data imputation
      tmp.gene <-  data.frame(apply(tmp.gene, 2, as.numeric))
      miss.idx <- which(colSums(tmp.gene ==0) >0)
      if(length(miss.idx) >0){
        for(i in colnames(tmp.gene)[ miss.idx ]) {
          tmp.index = which(tmp.gene[,i] == 0)
          if(length(tmp.index)>1){
            tmp.gene[tmp.index,i]  = mean(tmp.gene[,i], na.rm = T)
          }
        }
      }
      
      # -------------------------------------------- Marginally t-test and results 
      result.p <- c() ; result.t <- c()
      for(i in colnames(tmp.gene)) {
        snpvalue <- tmp.gene[,i]
        if(length(unique(snpvalue))>1){
          result <- t.test(snpvalue ~ Label.sub,   var.equal = TRUE)
          result.p <- c(result.p, result$p.value )
          result.t <- c(result.t, result$statistic)
        }else{
          result.p <- c(result.p, NA )  ; result.t <-c(result.t, NA )
        }
      }
      result.df <- as.data.frame(colnames(tmp.gene))
      result.df$pvalue <- result.p
      result.df$tvalue <- result.t
      
      setwd(dir.output)
      write.csv(result.df , paste0("WholeGenome_MAF", mafthre, "_chr", chr, "_split", 
                                   l.idx,"_pvalue.csv"))
    }
    
  }else{ # chr >=7 
    selectsnp <- bed.dim.list[[chr]][2]
    # -------------------------------------------- load in data
    setwd(dir.data)
    tmp <-  read.plink(bed=paste0("MAF", mafthre,"_CHR" ,chr,".bed"),
                       bim=paste0("MAF", mafthre,"_CHR", chr,".bim"),
                       fam=paste0("MAF", mafthre,"_CHR", chr,".fam"),
                       na.strings = c("0", "-9"), sep = "." ,
                       select.subjects = NULL, select.snps = selectsnp)
    tmp.gene <- as.data.frame(as.matrix(tmp$genotypes))  
    tmp.gene <- tmp.gene[order(match(rownames(tmp.gene),subject.ID)), ]
    rm(list=c("tmp", "tmp.gene"))
    
    # -------------------------------------------- Missing genetype data and imputation
    tmp.gene <-  data.frame(apply(tmp.gene, 2, as.numeric))
    miss.idx <- which(colSums(tmp.gene ==0) >0 )
    if(length(miss.idx) >0){
      for(i in colnames(tmp.gene)[ miss.idx ]) {
        tmp.index = which(tmp.gene[,i] == 0)
        if(length(tmp.index)>1){
          tmp.gene[tmp.index,i]  = mean(tmp.gene[,i], na.rm = T)
        }
      }
    }
    
    # -------------------------------------------- Marginally t-test and results 
    result.p <- c() ; result.t <- c()
    for(i in colnames(tmp.gene)) {
      snpvalue <- tmp.gene[,i]
      if(length(unique(snpvalue))>1){
        result <- t.test(snpvalue ~ Label.sub,   var.equal = TRUE)
        result.p <- c(result.p, result$p.value )
        result.t <- c(result.t, result$statistic)
      }else{
        result.p <- c(result.p, NA )  ; result.t <-c(result.t, NA )
      }
    }
    result.df <- as.data.frame(colnames(tmp.gene))
    result.df$pvalue <- result.p
    result.df$tvalue <- result.t
    
    setwd(dir.output)
    write.csv(result.df , paste0("WholeGenome_MAF", mafthre, "_chr", chr, "_pvalue.csv"))
  }
}

rm(list = ls(pattern="^tmp"))
rm(list = c("chr", "result.df", "miss.idx", "result.p", "snpvalue" ))




## 3.2 Combine all the marginal test results across different chrs -------
setwd(dir.output)
for (chr in 1:22){ print(chr)
  if(chr <7){
    tmp.size <- bed.dim.list[[chr]][2]
    l <- split(1:tmp.size, ceiling(seq_along(1:tmp.size)/100000))
    for (l.idx in 1:length(l)){
      tmp <- 
        read.csv(paste0("WholeGenome_MAF", mafthre, "_chr", chr, "_split", l.idx,"_pvalue.csv"))
      tmp$chr = chr
      if(l.idx == 1 & chr == 1){
        strong.test = tmp
      }else{ strong.test = rbind(strong.test, tmp)}
    }
  }else{
    tmp <- read.csv(paste0("WholeGenome_MAF", mafthre, "_chr", chr,"_pvalue.csv"))
    tmp$chr = chr
    strong.test = rbind(strong.test, tmp)
  }
}
rm(list = ls(pattern="^tmp"))
rm(list = c( "chr", "l" , "l.idx", "i", "result", "result.t"))
dim(strong.test) 


# 3.3 Select marginally strong SNPs based on the pre-determined p value threshold --- 
strong.test.select <-strong.test[strong.test$pvalue < pvalue.strong.thred, ]
dim(strong.test.select)[1]  
rm(list =c("strong.test"))

# Create a empty list to store SNPs within each chr. 
# Note: chr.snp is a list, each element of this list contains all the marginal strong SNPs in that chr 
chr.snp <- list() 


## 3.4 Add additional SNPs based on prior knowledge (if needed)  ----------
# Note: In our study, we add additional SNPs from the reference paper by Fritsche, Lars G., et al. "A large genome-wide association study of age-related macular degeneration highlights contributions of rare and common variants." Nature genetics 48.2 (2016): 134-143.

# For example, if you want to add the following SNPs on Chr1 and Chr3 into our marginally strong SNPs list
# For chr1, add the below SNPs
chr.snp[[1]]= c("rs10922109", "rs570618", "rs121913059", "rs148553336", "rs187328863", "rs61818925", "rs35292876", "rs191281603")
# For chr3, ad the below SNPs
chr.snp[[3]] =c("rs140647181", "rs55975637")


## 3.5 Combine the marginally strong SNPs together  ------- 
# marginally strong SNPs = marginally test detected SNPs + Prior-knowledge suggestes SNPs (if any)
a = 0 ;b = 0 ; c = 0
for( chr in 1:22){
  tmp <- chr.snp[[chr]]
  snp.ttest <- strong.test.select[strong.test.select$chr==chr, ]$colnames.temp.gene.
  a <- a + sum(tmp %in% as.character(snp.ttest))
  b <- b + length(as.character(snp.ttest))

  tmp <- unique(c(tmp, as.character(snp.ttest)))
  # Delete SNPs with strange names 
  chr.snp[[chr]] <- noquote(tmp[which(startsWith(tmp, "rs"))])

  c <- c + length(chr.snp[[chr]] )
}

# Summary information 
print(paste("Number of added strong SNPs overlap with the marginally strong SNPs by marginally tests is ", a))
print(paste("Number of marginally strong SNPs by marginally tests is ", b))
print(paste("Number of all marginally strong SNPs is ", c))

rm(list = ls(pattern="^tmp"))
rm(list = c("a", "b", "c", "snp.ttest"))


## 3.6  Summarize marginally strong SNP information -----
tmp.strong <- c()
tmp.chr <- 
for( chr in 1:22){
  tmp.strong <- c(tmp.strong, chr.snp[[chr]] )
  tmp.chr <- c(tmp.chr, rep(chr, length(chr.snp[[chr]])))
}
 
strong.info <- as.data.frame(cbind(  tmp.chr,  tmp.strong) )
colnames(strong.info) <- c("chr", "snp")
rownames(strong.info) <- NULL
strong.info$snp <- as.character(strong.info$snp)
head(strong.info)
rm(list = ls(pattern="^tmp"))
rm(list = c("chr"))




# 4. Search for marginally strong SNP correlated SNPs based on LD -------
# In this section, using the MAF-filtered data set, we first calculated the LD for each of the marginal strong SNP with other SNPs on the same chr.  R square measurement of LD is used. The SNPs pass the pre-determined LD threshold (LD.thre) will be returned . 
 
## 4.1  Calculate LD matrices for each marginal strong SNP with other SNPs on the same chr -------

LD.matrix.list <- list()  
result.n.snp.list <- list()  
result.n.strong.snp.list <- list()  
index.chr  <- c() # The snps within which chr. 

for(chr in which(lengths(chr.snp)!=0) ){
  # chr number 
  print(paste("chr", chr, chr.snp[[chr]])) 
  
  # How many marginal strong on this chr.
  print(paste("length", length(chr.snp[[chr]]))) 
  
  # Find the SNPs on the same chr with the marginally strong SNP, and with LD > LD.thred
  i= 1
  for(snpi in chr.snp[[chr]]){ print(snpi)
    setwd(dir.data)
    cmd <- str_c("plink --bfile MAF0.1_CHR",chr,"  --r2  --ld-window-r2 ", LD.thre,
                 " --ld-window 15000000 --ld-window-kb 20000000000",
                 " --ld-snp " ,snpi,
                 " --out WholeGenomePlus52snps_", LD.thre, "_LD_Chr",chr)
    system(cmd)
    try <- read.table(paste0("WholeGenomePlus52snps_", LD.thre, "_LD_Chr",chr,".ld"), 
                      sep="", h = T, as= T)
    
    # --------------------------------------------------- LD matrix format 
    tmp.matrix <- matrix(0, nrow = length(unique(c(try$SNP_A, try$SNP_B)))
                         , ncol = length(unique(c(try$SNP_A, try$SNP_B))))
    row.names(tmp.matrix) = unique(c(try$SNP_A, try$SNP_B))
    colnames(tmp.matrix) = unique(c(try$SNP_A, try$SNP_B))
    if (dim(try)[1] > 0 ){
      for(i.dim in 1:dim(try)[1]){
        i.row = which(unique(c(try$SNP_A, try$SNP_B)) == try$SNP_A[i.dim])
        i.col = which(unique(c(try$SNP_A, try$SNP_B)) == try$SNP_B[i.dim])
        tmp.matrix[i.row, i.col] =  try$R2[i.dim]
        tmp.matrix[i.col, i.row] =  try$R2[i.dim]
      }
    } 
    LD.matrix.list[[paste0("chr", chr, "snp", snpi)]] <- tmp.matrix
    if(i==1){
      try.combine = try
      i = 0
    }else{
      try.combine = rbind(try.combine, try)
    } 
  }
  
  
  # --------------------------------------------------- summary 
  # Number of SNPs detected 
  result.n.snp.list[[chr]] <-length(unique(c(try.combine$SNP_A, try.combine$SNP_B)))
  # Number of marginally strong SNPs with LD-correlated SNPs (LD > LD.thred)
  result.n.strong.snp.list[[chr]] <- sum(unique( c(try.combine$SNP_A, 
                                                 try.combine$SNP_B)) %in% chr.snp[[chr]])
  
  # --------------------------------------------------- index.chr 
  if(result.dim.list[[chr]]>0){  index.chr<-  c(index.chr , chr)}
}


rm(list = ls(pattern="^tmp"))
rm(list = ls(pattern="^i"))
rm(list=c("try", "try.combine", "chr", "snpi", "cmd" ))



## 4.2 Summarize the LD-correlated SNPs------- 
LDcorrelated.snp.info <- 
  data.frame("chr" = index.chr, 
             # Number of SNPs detected in each chr 
             "n.snp" = Reduce(rbind, result.n.snp.list) ,  
             # Number of marginally strong SNPs with LD-correlated SNPs  in each chr 
             "n.strong.snp" = Reduce(rbind, result.n.strong.snp.list))
rownames(LDcorrelated.snp.info) <- NULL 
LDcorrelated.snp.info
rm(list = c("result.dim.list", "result.unique.snp.list"))

##  LD.matrix.list is a list, each element contains the LD matrix for a marginally strong SNPs with all SNPs on the same chr (Pass MAF and LD thresholds). 


## 4.3  Summarize marginally strong SNP information (Result is the same as 3.6)  
# strong.info <- names(LD.matrix.list)
# strong.info <- as.data.frame(Reduce(rbind, strsplit(strong.info, "snp"))) 
# colnames(strong.info) <- c("chr", "snp")
# rownames(strong.info) <- NULL
# strong.info$snp <- as.character(strong.info$snp)
# head(strong.info)
# # chr         snp
# # 1 chr1  rs10922109
# # 2 chr1    rs570618
# # 3 chr1 rs121913059
# # 4 chr1 rs148553336
# # 5 chr1 rs187328863
# # 6 chr1  rs61818925



## 4.3 Further filter correlated SNPs with LD > alpha ------------------------- 
# Note: We can further filter the LD-correlated SNPs using a pre-determined threshold (alpha, which is larger than the LD.thred). The reason why we include both section LD.thred and alpha is because users might tune different alpha level to get better prediction results.  
LD.matrix.list.copy <- LD.matrix.list
for(i in 1:length(LD.matrix.list)){
  tmp.matrix <- LD.matrix.list[[i]]
  tmp.matrix[ abs(tmp.matrix) < alpha ] <- 0L
  tmp.matrix <- tmp.matrix[which(colSums(tmp.matrix) != 0) , which(colSums(tmp.matrix) != 0)]
  LD.matrix.list[[names(LD.matrix.list)[i]]] <- tmp.matrix
}
rm(list=c("tmp.matrix", "i"))
# names(LD.matrix.list) # e.g. "chr6snprs9273481" 




# 5. Construct SNP clusters based on LD matrices----------------------------
# In this section, SNP clusters are constructed, each of which contains at least one marginally strong SNP and all the LD-correlated SNPs. We combine the clusters if they have one or more SNPs in common. 

## 5.1 Search SNP clusters ---------------------------------------------- 
# To make selected.snps.name $chr$strong signal  --> all snps in each cluster of strong signals
cluster.snp.list <- list()

# Count how many strong signals do and do not have LD.correlated SNPs 
tmp.count.no = 0 ; # strong signal have no cluster
tmp.count.yes = 0 ; # strong signal have cluster

chr = "initiation" 
tmp.chr.list <- list()
for (tmp.strongsnp in names(LD.matrix.list)){
  tmp.LDmatrix  <- LD.matrix.list[[tmp.strongsnp]]
  tmp.chr <-  strsplit(tmp.strongsnp, "snp") [[1]][1]
  tmp.strongsnp <- strsplit(tmp.strongsnp, "snp") [[1]][2]
  
  if (chr!= tmp.chr) { tmp.chr.list <- list()}
  
  if(is.null(dim(tmp.LDmatrix))){  
    tmp.chr.list[[tmp.strongsnp]] <- tmp.strongsnp
    tmp.count.no = tmp.count.no + 1
  }else{
    tmp.name <- colnames( tmp.LDmatrix )
    tmp.chr.list[[tmp.strongsnp]] <- tmp.name[which(startsWith(tmp.name , "rs"))]
    tmp.count.yes = tmp.count.yes + 1
  }
  cluster.snp.list[[tmp.chr]]  <- tmp.chr.list
  chr = tmp.chr 
}

cluster.snp.list <- cluster.snp.list[lengths(cluster.snp.list ) > 0]

print(paste("N of marginally strong SNP have associated cluser (at least one LD-correlated SNP) is ", tmp.count.yes))
print(paste("N of marginally strong SNP does not have asscociated cluser is ", tmp.count.no))
rm(list = ls(pattern="^tmp"))




## 5.2  Combine the SNP clusters into one if they share commone SNP(s) ----
# Function for this aim 
unique.component <- function(componentlist  ){
  componentlist <-  componentlist[lengths(componentlist) > 0L]
  if(length(unlist(componentlist)) == length (unique( unlist(componentlist)))){
    component  <- c() ; group <- c()
    for(i in 1:length( componentlist)){
      component_current <- componentlist[[i]]
      group_current <- rep(i, length(component_current))
      component <- c(component , component_current) ;
      group <- c(group, group_current)
    }
    component.group.df <- as.data.frame(cbind(component = component, group= group))
    return(component.group.df)
  }else{
    B = length(componentlist)
    for(a in 1:(B -1)){
      avec <- componentlist[[a]]
      if ( length(avec) >0 ){
        for( b in  (a+1):B ) {
          bvec <- componentlist[[b]]
          if(length(intersect(avec, bvec)) !=0 ){
            avec = unique(c(avec, bvec));
            componentlist[[b]] <- character(0)
          }
        }
      }
      componentlist[[a]] <- avec
    }
    componentlist <-  componentlist[lengths(componentlist) > 0L]
    unique.component(componentlist)
  }
}


final.cluster.list <- list()
final.cluster.dimbychr <- list(); # list format 
final.cluster.dimbychr2 <- list() # df format 
for( chr in names(cluster.snp.list) ){ print(chr)
  tmp.thischr <- cluster.snp.list[[chr]]
  if( !is.null(tmp.thischr)  & length(cluster.snp.list[[chr]]) >0 ){     
    # combine clusters if have common snps
    final.cluster.list[[chr]] <-  unique.component(componentlist = cluster.snp.list[[chr]] )
    # the unique cluster of each rs
    tmp.chr.group = unique( final.cluster.list[[chr]][,"group"])
    # N of rs in each cluster
    tmp.sum <- c()
    for( tmp.group in tmp.chr.group){
      tmp.sum <- c(tmp.sum, length(which( final.cluster.list[[chr]][,"group"] == tmp.group)))
    }
    final.cluster.dimbychr[[chr]] <- tmp.sum
    final.cluster.dimbychr2[[chr]] <- cbind( rep( chr, length(tmp.sum)), tmp.sum)
  }
}
rm(list = ls(pattern="^tmp"))
rm(list = c("chr"))



## Below are examples for what 'final.cluster.dimbychr' and 'final.cluster.dimbychr2' looks like 
# final.cluster.dimbychr
# $chr1
# [1]  22 351 102   1   1  34  87   2  49  64 101   2   1   1  10   3   7   2   1
# [20]   1  28   1   1  34   1   1   1   2   1   2   2   1   2   1  38   3   2   1
# [39]  22   3 251   2   1   1  58  30   1  22   1   2  12   1   1   1 110   2   2
# [58]   1   2 133   1   1   1   1   1  61   2   4  16

# final.cluster.dimbychr2
# [86,] "chr6" "288"  
# [87,] "chr6" "2"    
# [88,] "chr6" "1"    
# [89,] "chr6" "1"    
# [90,] "chr6" "1"    
# [91,] "chr6" "1"    
# [92,] "chr6" "26"   



## 5.3 Information on final SNP clusters-------------------
### Descriptives -----
# N of cluster contains only one SNP (1 marginally strong SNP)
sum(unlist(final.cluster.dimbychr) ==1 ) 
# N of cluster have more than one SNP
sum(unlist(final.cluster.dimbychr) > 1 ) 
# sum(unlist(final.cluster.dimbychr)[unlist(final.cluster.dimbychr) >1 ])  
# sum((unlist(final.cluster.dimbychr)))  
rm(final.cluster.dimbychr)

# Summarize the size of each cluster -- using final.cluster.dimbychr2
final.cluster.dimbychr2 <-as.data.frame( Reduce(rbind, final.cluster.dimbychr2 )) 
names(final.cluster.dimbychr2) <- c("chr", "size")
final.cluster.dimbychr2$size <- as.numeric(as.character(final.cluster.dimbychr2$size))

dim(final.cluster.dimbychr2) # 410  2
final.cluster.dimbychr2$chr = as.factor(final.cluster.dimbychr2$chr )
final.cluster.dimbychr2$chr   <- ordered(final.cluster.dimbychr2$chr ,
                                         levels = c("chr1" ,"chr2"  , "chr3" , "chr4" , 
                                                    "chr5" , "chr6" , "chr7",  "chr8", 
                                                    "chr9", "chr10" ,"chr12" ,"chr13",
                                                    "chr14", "chr15", "chr16" ,"chr17" ,
                                                    "chr19","chr20", "chr22" )
)
# N of SNPs within all clusters 
sum(final.cluster.dimbychr2$size)  
# Five number of size of all clusters
fivenum(final.cluster.dimbychr2$size)  
# Max size of all clusters
final.cluster.dimbychr2$size[which.max(final.cluster.dimbychr2$size) ]  


### Visulization -----
# library(ggplot2)
# Histrogam of cluster sizes by chr
p <-  final.cluster.dimbychr2  %>%
  filter( final.cluster.dimbychr2$size > 1) %>%
  # mutate(text = fct_reorder(chr, value)) %>%
  ggplot( aes(x=size, color=chr, fill=size)) +
  geom_histogram(alpha=0.6, binwidth = 1) +
  ggtitle(paste("Histrogam of cluster sizes by chr \n LD =", alpha)) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Size") +
  facet_wrap(~chr)
p

# Histogram of cluster sizes in all chr
hist( final.cluster.dimbychr2$size, xlab="",
      main = paste("Histrogam of cluster sizes in all chr \n LD =", alpha) )





# 6 Prepare data which contains all SNPs in the detected SNP clusters  ------------- 
# In this section, we will prepare our final analysis data set, which contains all the SNPs we detected in the clusters. Missing genotype data imputation is performed when needed.

## 6.1 Extract data -------  
data <- NULL
set(dir.data)
for(chr in 1:22){  print(chr)
  tmp <- final.cluster.list[[paste0("chr", chr)]]
  tmp$group<- as.numeric(as.character(tmp$group))
  tmp.snpnames <- as.character(tmp$component) 
  
  # Note: Similar to Section 2, due to computation time and memory concern, I split the larger chr1-6 into smaller blocks. 
  if(chr < 7 ){
    tmp.size <- dim.list[[chr]][2]
    l <- split(1:tmp.size, ceiling(seq_along(1:tmp.size)/100000))
    for(l.idx in 1:length(l)){
      selectsnp <- l[[l.idx]]
      # -------------------------------------------- load in data
      tmp <-  read.plink(bed=paste0("MAF", mafthre,"_CHR" ,chr,".bed"),
                         bim=paste0("MAF", mafthre,"_CHR", chr,".bim"),
                         fam=paste0("MAF", mafthre,"_CHR", chr,".fam"),
                         na.strings = c("0", "-9"), sep = "." ,
                         select.subjects = NULL, select.snps = selectsnp)
      tmp.gene <- as.data.frame(as.matrix(tmp$genotypes)) # head(  tmp.gene)
      tmp.rownames <- row.names(tmp.gene)
      tmp.gene <- tmp.gene[which(tmp.rownames %in% subject.ID), ]
      tmp.gene <- tmp.gene[order(match(rownames(tmp.gene),subject.ID)), ]
      rm(list=c("tmp", "tmp.rownames"))
      
      # --------------------------------------------  snps
      tmp.snpname <-  tmp.snpnames[which(tmp.snpnames %in% colnames(tmp.gene))]
      if(length(tmp.snpname ) >=1){
        tmp.gene <- as.data.frame(tmp.gene[, tmp.snpname  ])
        colnames(tmp.gene) <- tmp.snpname
        # combine
        if(is.null(data )  ){
          data  <-  tmp.gene
        }else    {
          data = cbind(data, tmp.gene)
        }
      }
      
    }
    
  }else{ # chr >=7
    # -------------------------------------------- load in data
    tmp <-  read.plink(bed=paste0("MAF", mafthre,"_CHR" ,chr,".bed"),
                       bim=paste0("MAF", mafthre,"_CHR", chr,".bim"),
                       fam=paste0("MAF", mafthre,"_CHR", chr,".fam"),
                       na.strings = c("0", "-9"), sep = "." ,
                       select.subjects = NULL, select.snps = NULL)
    tmp.gene <- as.data.frame(as.matrix(tmp$genotypes)) # head(  tmp.gene)
    tmp.rownames <- row.names(tmp.gene)
    tmp.gene <- tmp.gene [which(tmp.rownames %in% subject.ID), ]
    tmp.gene <- tmp.gene [order(match(rownames(tmp.gene),subject.ID)), ]
    rm(list=c("tmp", "tmp.rownames"))
    # --------------------------------------------  snps
    tmp.snpname <-  tmp.snpnames[which( tmp.snpnames  %in% colnames( tmp.gene))]
    if(length(tmp.snpname) >=1){
      data = cbind(data, tmp.gene[, tmp.snpname])
    }
  }
  
}
rm(list = ls(pattern="^tmp"))
rm(list = c( "l", "l.idx"))


## 6.2 Missing genotype data imputation  ------------------------------------------------------- 
for(i in 1:dim(data)[2]){
  data[,i] <- as.numeric(as.character(data[,i]))
  if(sum(is.na(data[,i]) >0 )) {
    data[which(is.na(data[,i])), i] <- mean(data[,i], na.rm=T)
  }
}
rm(i)



## 6.3 Check  ----- 
# Delete SNPs without "rs" number if needed. 
which(startsWith( colnames(data) , "rs")==F)  





# 7. Detect marginally weak signals based on cluster-adjusted effect sizes ---------------
# In this section, we calculate the SNP cluster-adjusted effect size (iff condition value as in mLDA) for each SNP with in each cluster. Details please refer to Li, Y., Hong, H.G. and Li, Y., 2019. Multiclass linear discriminant analysis with ultrahigh‐dimensional features. Biometrics, 75(4), pp.1086-1097. Then, based on this value we selected the marginally weak SNPs set. The final informative SNP sets contains all the marginally weak and strong SNPs.



##  7.1 cluster-adjusted effect sizes (iff condition value as in mLDA)  ------------- 

# Function for calculating the precision matrix, or the inverse of covariance matrix 
# library(matrixcalc)
my.inv <- function(X, eps=1e-12){
  if (is.singular.matrix(X) ){
    eig.X <- eigen(X, symmetric=TRUE)
    P <- eig.X[[2]]
    lambda <- eig.X[[1]]
    ind <- lambda > eps
    lambda[ind] <- 1/lambda[ind]
    lambda[!ind] <- 0
    ans <- P%*%diag(lambda, nrow=length(lambda))%*%t(P)
    return(ans)
  }else{
    return(solve(X))
  }
}



# cluster-adjusted effect sizes
iffcond.list <- list() # a list. Each element of the list is info for one marginally strong SNP. Element name consists information of the chr and the Cluster/Group the SNP is located. Element contents are for the associated cluster. Information includes, the rank of cluster-adjusted effect size (in absolute value) of the marginally strong SNP within the cluster (decreasing order),  chr number, total n. of snps in the cluster, n of marginal strong SNP


iffcond.df  <-  list() # a list. Each element of the list is info for one marginally strong SNP Element name consists information of the chr and the Cluster/Group the SNP is located.  Element content is a data.frame for all the SNPs within the associated cluster. Information includes, chr number, group number, SNP rsnumber, marginally weak or strong SNP, cluster-adjusted effect size(iff condition value) ,n of snp in this cluster. 

sum.chr <-c() ;   
sum.ncluster <- c()
sum.nifflarge <- c()
sum.strong <- c()


for(chrname in names(final.cluster.list)){ print(chrname)
  chr =  as.numeric(as.character(Reduce(rbind, strsplit(chrname, "chr"))[2]))
  
  tmp <- final.cluster.list[[chrname]]
  tmp$group <- as.numeric(as.character(tmp$group))
  chr.B = unique(tmp$group)
   
  for(tmp.chr.B.index in chr.B){
    tmp.snpnames <-  as.character(tmp[tmp$group == tmp.chr.B.index, 1] )
    
    # For cluster only contains one SNP
    if(length(tmp.snpnames) == 1){
      
          # Summary version 1
          tmp.combine.update <- c( "NULL", # the rank of cluster-adjusted effect size (in absolute value) of the marginally strong SNP within the cluster (decreasing order). Because there is only 1 strong SNP, we does not care. 
                                   chr,    # chr
                                   1 ,     # total n of snps in the cluster
                                   1       # n of strong
          )
          names(tmp.combine.update) <- tmp.snpnames
          iffcond.list[[paste("chr", chr, "group", tmp.chr.B.index) ]] <- tmp.combine.update
            
          # Below is for summary version 2, used later 
            sum.chr <- c(sum.chr, chr)             # chr.
            sum.ncluster <- c(sum.ncluster, 1)     # n of SNPs in the B
            sum.nifflarge <- c(sum.nifflarge, NA ) # n of weak SNPs with iff cond value large than that of strong. NA means no weak SNP in the cluster 
            sum.strong <- c(sum.strong, 1 ) # n of strong 
            
          # Below is for summary version 3,
            # Calculate iff condition values 
            tmp.gene <- data[, tmp.snpnames] 
            tmp.meanDiff.s <- mean(tmp.gene[which( Label.sub==1)  ]) -  mean(tmp.gene[which( Label.sub==1)  ])
            tmp.iff <- as.vector(1/(var(tmp.gene )) * tmp.meanDiff.s)
            iffcond.df[[paste("chr", chr, "group", tmp.chr.B.index)]] <-  
              # chr, chr_group,       snp #rs,       s/w,     iff condition value, n.snp
              c(chr, tmp.chr.B.index, tmp.snpnames, "strong" , tmp.iff, 1 ) 
    
      
    # For clusters contain more than one SNPs
    }else{  
      # Calculate iff condition values 
      tmp.gene <- data[,tmp.snpnames]
      tmp.colnames  <- colnames( tmp.gene)
      tmp.meanDiff.s <- apply(as.matrix(tmp.gene[which( Label.sub==1),]), 2, mean) - apply(as.matrix(tmp.gene[which(Label.sub==2),]), 2, mean)
      tmp.combine <- as.vector(my.inv(cov(tmp.gene )) %*%  tmp.meanDiff.s)
      names(tmp.combine ) = tmp.colnames
      # Compare all SNPs to marginal strong iff
      ranks <- order(-abs(tmp.combine))
      tmp.combine <- tmp.combine[ranks]
      tmp.n <- which(names(tmp.combine) %in% marginalstrong)
      
      # Summary version 1
      iffcond.list[[paste("chr", chr, "group", tmp.chr.B.index) ]] <- 
        c(tmp.combine,
          chr,
          length(ranks ) , # total snps in the B
          tmp.n            # of strong
        )
      
      # Below is for summary version 2, used later 
      sum.chr <- c(sum.chr, chr)
      sum.nifflarge.tmp <-tail(tmp.n, n=1) - length(tmp.n)
      sum.ncluster <- c(sum.ncluster, length(ranks )) # total snps in the B
      sum.nifflarge <- c(sum.nifflarge,  sum.nifflarge.tmp ) #  of strong
      sum.strong <- c(sum.strong, length(tmp.n) )
      
      # Summary version 3
      tmp.sw <- rep("weak", length(ranks))
      tmp.sw[tmp.n] <- "strong"
      
      tmp.df <- cbind( rep(chr, length(ranks )), 
                       rep(tmp.chr.B.index,length(ranks ))  ,
                       names(tmp.combine[1:length(ranks )]),
                       tmp.sw, 
                       tmp.combine[1:length(ranks )], 
                       rep(length(ranks ), length(ranks )))
      iffcond.df[[paste("chr", chr, "group", tmp.chr.B.index)]] <-  tmp.df
      
    }
  }
}

rm(list = ls(pattern="^tmp"))



# Summarize the information for each  (version 2)
iffcond.df <- as.data.frame(Reduce(rbind, iffcond.df))
rownames(iffcond.df) <- NULL
colnames(iffcond.df) <- c("chr", "B", "snpname", "sw", "iff", "ncluster")
iffcond.df$chr <- as.numeric(as.character(iffcond.df$chr))
iffcond.df$B <- as.numeric( as.character(iffcond.df$B))
iffcond.df$sw <- as.character(iffcond.df$sw) # Marginally weak or strong SNP
iffcond.df$iff <- as.numeric( as.character( iffcond.df$iff))
iffcond.df$ncluster <- as.numeric(as.character( iffcond.df$ncluster ))
#   chr B    snpname     sw  iff ncluster
# 1   1 1 rs10922106 strong  0.79       27
# 2   1 1 rs10922105 strong .3740       27
# 3   1 1 rs10922108 strong .7677       27
# 4   1 1 rs10737680 strong .3624       27
# 5   1 1  rs6688272 strong .7618       27
# 6   1 1  rs1410996 strong .3612       27
iffcond.df$absiff <- abs(iffcond.df$iff) 

# Number of marginally weak and strong SNPs
table(iffcond.df$sw)
# Number of group in each chr
tapply(iffcond.df$B, iffcond.df$chr, max)
# Summary of iff condition  
iff.sum <- as.data.frame(sum.chr)
# N of SNPs within each cluster
iff.sum$sum.ncluster <-  sum.ncluster   
# N of clusters with weak signal iff larger than strong signal in the same cluster. NA means no weak signals in this cluster.
iff.sum$sum.nifflarge <- sum.nifflarge  
iff.sum$sum.strong <- sum.strong

iff.sum






## 7.2 Select final SNPs  ------------------------------------
# In this section, we select the informative SNPs, which contains the marginally strong SNPs (as in Section 3) and marginally weak SNPs based on the cluster-adjusted effect sizes in Section 6. 

## Select marginally weak informative SNPs
# Note: For example, we can select the top 1000 weak signals which have the largest cluster-adjustd effect sizes (in absolute value)
select.weak.iff <- iffcond.df %>% 
  filter(sw =="weak") %>% 
  arrange(desc(absiff))  %>% 
  slice(1:1000)

# chr and B (group indicator) are the marginally weaks SNPs locates 
weak_B_num <- select.weak.iff %>% select(chr, B) %>% 
  group_by(chr) %>% mutate(nB = length(unique(B)))  %>%
  select(chr, nB) %>% as.data.frame
weak_B_num <- weak_B_num[duplicated(weak_B_num)==F, ]  
weak_B_num <-weak_B_num[order(weak_B_num$chr), ]
plyr::tapply(select.weak.iff$B, select.weak.iff$chr, count)


## Select all marginally strong informative SNPs
select.strong.iff<- iffcond.df %>% 
  filter(sw =="strong") %>%
  group_by(chr, B) %>% 
  slice(n=1)


###  Final select marginally strong and weak SNPs
select <- rbind(select.strong.iff, select.weak.iff)





# 8 Outcome classification based on selected SNPs  ---------------------------
# In this section, we build outcome classification model based on 8.1 MLDA decision rule or Fisher's discriminant decision rule; 8.2 Ridge regression. To compare the effect of weak signals on outcome prediction, prediction model is built using   1) all the marginally strong SNPs, and 2) both the marginally strong and weak SNPs. More specifically, we split the corresponding data set into 5 folds. Using 4 folds to build the model and predict on the left-fold. Repeat this process 5 times so that every data point get a prediction. The performance of prediction is evaluated using area under the ROC Curve. 


library(caret)
library(glmnet)
library(PRROC)
library(pROC)
 
dim(data)

# all strong snp names
strong.snp <- iffcond.df.orignal[sw=="strong",]$snpname
# all selected weak and strong snp names
strongNweak.snp<- select$snpname

# iffcond.df.orignal <- iffcond.df




## 8.1 Build model and Prediction based on mLDA prediction rule---- ---------------------------------
# Note: descision rule is based on Fisher's discriminant rule 



seed.i = 123 ;  nFolds = 5
set.seed(seed.i)

# function for cross-validation 
summary_mLDArule <- function(nFolds = 5 , X, y, Z  ){
  
  folds.l1  <- createFolds(y, k = nFolds)  # outer is level1 
  pred.y.mlda <- c(); truth.y  <- c()
  
  
  Omega.cis.mat =my.inv(cov(as.matrix(X )))
  ZOmega <- qr.solve(cov(Z))
  
  folds.l1  <- createFolds(y, k = nFolds)  # outer is level1 
  pred.y.mlda <- c(); truth.y  <- c()
  for(i.l1 in 1:nFolds){ #i.l1 = 1 
    print(i.l1)
    
    ##################################### split data 
    index.real.l1out =folds.l1[[i.l1]]
    index.real.l1in = setdiff(seq_len(length(y)),  index.real.l1out)
    
    X_Train = X[index.real.l1in, ] ; y_Train = y[index.real.l1in]
    X_Test = X[index.real.l1out, ]; y_Test = y[index.real.l1out]
    
    Z_Train = Z[index.real.l1in, ]  
    Z_Test = Z[index.real.l1out, ] 
    
    X.new.s = cbind(Z_Test, X_Test)
    
    ########################################## Genotype data 
    Omega.cis.mat =  Omega.cis.mat
    # Omega.cis.mat =my.inv(cov(as.matrix(X_Train)))
    meanDiff = apply(as.matrix(X_Train[which(y_Train==2),]), 2, mean) - apply(as.matrix(X_Train[which(y_Train==1),]), 2, mean)
    meanAvg <- (apply(as.matrix(X_Train[which(y_Train==2),]), 2, mean) + apply(as.matrix(X_Train[which(y_Train==1),]), 2, mean))/2
    pXs2 <- length(meanDiff)
    # Fisher <- as.matrix((X_Test - mean.cis.mat/2))%*%Omega.cis.mat%*% mean.cis.mat
    
    ################################ Include covariate information  Z
    ZOmega =  ZOmega
    ZmeanDiff <- apply(as.matrix(Z_Train[which(y_Train==2), ]), 2, mean) - apply(as.matrix(Z_Train[which(y_Train==1), ]), 2, mean)
    ZmeanAvg <-( apply(as.matrix(Z_Train[which(y_Train==2), ]), 2, mean) + apply(as.matrix(Z_Train[which(y_Train==1), ]), 2, mean))/2
    # ZOmega <- qr.solve(cov(Z_Train))
    pZ <- dim(Z_Train)[2]
    
    ################################# combine phenotype and covariate info together
    meanAvg.s <- c(ZmeanDiff,  meanAvg)
    meanDiff.s2 <- c(ZmeanDiff,  meanDiff)
    Omega.s2 <- matrix(0, (pZ+pXs2), (pZ+pXs2))
    Omega.s2[1:pZ, 1:pZ] <- ZOmega
    Omega.s2[(pZ+1):(pZ+pXs2), (pZ+1):(pZ+pXs2)] <- Omega.cis.mat
    
    Fisher <-   as.matrix(X.new.s-meanAvg.s)%*%Omega.s2%*%meanDiff.s2
    PredClass <- (Fisher>=0)
    
    ################################## Prediction  Results 
    pred.y.mlda <- c(pred.y.mlda, PredClass)
    truth.y  <- c( truth.y , y_Test)
    
  }
  
  return( list(truth.y=truth.y, pred.y.mlda=pred.y.mlda) )
}


###  8.1.1 Prediction based on marginally strong SNPs alone-----

# Prepate data 
y = Label # outcome information 
X = data[, which(colnames(data) %in% strong.snp)]  # SNP information 
Z =  as.data.frame(cbind(Demog$age, Demog$sex)) ; # Covariates 
colnames(Z) = c("age", "sex") 

# Model 
set.seed(seed.i)
result_strong_together <-summary_mLDArule(nFolds = 5 , X, y, Z  )
truth.y = result_strong_together[["truth.y"]] - 1
pred.y.mlda= result_strong_together[["pred.y.mlda"]]

# Predictoin performance
mean(pred.y.mlda==  truth.y)   
auc(pred.y.mlda,  truth.y)  


## 8.1.2Prediction based on marginally strong and weak SNPs  ------------
# Prepare data 
y = Label  # outcome information 
X = data[,  which(colnames(data) %in% strongNweak.snp)]  # SNP information 
Z =  as.data.frame(cbind(Demog$age, Demog$sex))  # Covariates 
colnames(Z) = c("age", "sex")  
# Model 
set.seed(seed.i)
result_strong_together <-summary_mLDArule(nFolds = 5 ,X=X, y=y, Z=Z  )
truth.y = result_strong_together[["truth.y"]]-1
pred.y.mlda= result_strong_together[["pred.y.mlda"]]
# Prediction performance
mean(pred.y.mlda==  truth.y)  
auc(pred.y.mlda,  truth.y)  


#### 8.2 Build model and Prediction based on ridge regression --

##### 8.2.1 Prediction based on marginally strong SNPs alone-----
# Prepate data 
y = Label # outcome information 
X = data[, which(colnames(data) %in% strong.snp)]  # SNP information 

# model 
seed.i = 123 ;nFolds = 5
set.seed(seed.i)

folds.l1  <- createFolds(y, k = nFolds)  # outer is level1

pred.y.final.lasso <- c();
pred.y.final.ridge <- c(); truth.y  <- c()

for(i.l1 in 1:nFolds){ print(i.l1) #i.l1 = 1
  index.real.l1out =folds.l1[[i.l1]]
  index.real.l1in = setdiff(seq_len(length(y)),  index.real.l1out)
  
  
  X_Train = X[index.real.l1in, ] ; y_Train = y[index.real.l1in]
  X_Test = X[index.real.l1out, ]; y_Test = y[index.real.l1out]
  
  # ------------- ridge 
  print("ridge")
  model.ridge <-  glmnet(x= data.matrix(X_Train),   y= y_Train,  alpha = 0, family="binomial" , type.measure = "class", lambda = 1e-8)
  pred.pr <- predict(object= model.ridge, newx = data.matrix(X_Test), s =  1e-8,
                     type="response")
  pred.y <- as.numeric(pred.pr>0.5)
  
  pred.y.final.ridge <- c(pred.y.final.ridge, pred.y)
  truth.y  <- c( truth.y , y_Test)
  
  
}

truth.y = truth.y -1
# evaluation 
mean(pred.y.final.ridge==  truth.y)
auc(pred.y.final.ridge,  truth.y)


##### 8.2.2 Prediction based on marginally strong and weak SNPs  ----
# Prepare data 
y = Label  # outcome information 
X = data[,  which(colnames(data) %in% strongNweak.snp)] 

 

# model 
seed.i = 123 ;nFolds = 5
set.seed(seed.i)

folds.l1  <- createFolds(y, k = nFolds)  # outer is level1

pred.y.final.lasso <- c();
pred.y.final.ridge <- c(); truth.y  <- c()

for(i.l1 in 1:nFolds){ print(i.l1) #i.l1 = 1
  index.real.l1out =folds.l1[[i.l1]]
  index.real.l1in = setdiff(seq_len(length(y)),  index.real.l1out)
  
  
  X_Train = X[index.real.l1in, ] ; y_Train = y[index.real.l1in]
  X_Test = X[index.real.l1out, ]; y_Test = y[index.real.l1out]
  
  # ------------- ridge 
  print("ridge")
  model.ridge <-  glmnet(x= data.matrix(X_Train),   y= y_Train,  alpha = 0, family="binomial" , type.measure = "class", lambda = 1e-8)
  pred.pr <- predict(object= model.ridge, newx = data.matrix(X_Test), s =  1e-8,
                     type="response")
  pred.y <- as.numeric(pred.pr>0.5)
  
  pred.y.final.ridge <- c(pred.y.final.ridge, pred.y)
  truth.y  <- c( truth.y , y_Test)
  
  
}

truth.y = truth.y -1

# evaluation 
mean(pred.y.final.ridge==  truth.y)
auc(pred.y.final.ridge,  truth.y)






# 9 Permutation p value for each individual SNP ---------------------
# In this section, we use permutation test to get p-value for each SNP. Bascially, we permute the outcome label to preserve the cluster structure. We construct the null distribution for the cluster-adjusted effect size (iff-cond value) and calculate the probability of getting a value more extrame than observed cluster-adjusted effect size (iff-cond value). This process can be done through parallel computing. 


## 9.1 Function for permutation   ------ 
permute.label <- function(label, n , data, seed =111 ){ # n = Required number of permutations for a permutation-based p-value
  set.seed(seed)
  snpnames = colnames(data)
  
  # For cluster contains only 1 SNP
  if(length(snpnames)==1){
    meanDiff.s <- apply(as.matrix(data[which( label ==1),]), 2, mean) - 
      apply(as.matrix(data[which(label==2),]), 2, mean)
    omega <- 1/(var(data))
    original = as.vector(omega  *  meanDiff.s)
    names(original) = snpnames 
    for( i in 1:n){
      set.seed(seed+1+i)
      label.i <-  sample(label, length(label), FALSE)
      meanDiff.s <- apply(as.matrix(data[which( label.i ==1),]), 2, mean) - apply(as.matrix(data[which(label.i ==2),]), 2, mean)
      result.i <- as.vector( omega  * meanDiff.s)
      if(i==1){result = result.i}else{
        result = rbind(result, result.i)
      }
    }
    colnames(result) = snpnames
    
    result.p = c()
    for(snpi.in in snpnames){
      result.p = c(result.p ,  sum( abs(original[snpi.in]) <= abs(result[,snpi.in]) ) /(n)) 
    }
    names(result.p) <- snpnames
    
  # For cluster contains more than one SNPs
  }else{
    meanDiff.s <- apply(as.matrix(data[which( label ==1),]), 2, mean) - 
      apply(as.matrix(data[which(label==2),]), 2, mean)
    
    cov.df <- var(data)
    omega <- my.inv(cov.df )
    
    
    original = as.vector(omega  %*%  meanDiff.s)
    names(original) = snpnames 
    for( i in 1:n){
      set.seed(seed+1+i)
      label.i <-  sample(label, length(label), FALSE)
      meanDiff.s <- apply(as.matrix(data[which( label.i ==1),]), 2, mean) - 
        apply(as.matrix(data[which(label.i ==2),]), 2, mean)
      result.i <- as.vector( omega   %*%  meanDiff.s)
      if(i==1){result = result.i}else{
        result = rbind(result, result.i)
      }
    }
    colnames(result) = snpnames
    
    result.p = c()
    for(snpi.in in snpnames){
      result.p = c(result.p ,  
                   sum( abs(original[snpi.in]) <= abs(result[,snpi.in]) ) /(n)) # sum(abs(result[,snpi.in]) >= abs(original[snpi.in]))
    }
    names(result.p) <- snpnames
  }
  result <- result.p 
  return(result)
  
}

## 9.2 Permutation using parallel computing ------- 
library(foreach)
library(doParallel)
B = length(iffcond.list  )
cores= 32
cl <- makeCluster(cores)  # not to overload your computer
registerDoParallel(cl) 
single.p.list <- list()
single.p.list <- foreach(i.index= 1:B  )%dopar% {  # , .combine=cbind) %dopar% {
  
  tmp <- iffcond.list[i.index]
  tmp.listname <- names(tmp)
  tmp.snpname <- setdiff(names(tmp[[1]]), c("", NA))
  
  library(snpStats)
  library(glmperm)
  library(snpStats)
  library(matrixcalc)
  
  # -------------------------------------------- mLDA
  
  p.list.tmp <- list()
  tmp.data = data[, which(colnames(data) %in% tmp.snpname)]; 
  tmp.data =  tmp.data[, duplicated(t(tmp.data)) == FALSE]
  tmp.snpname1 =  tmp.snpname[duplicated(t(tmp.data)) == FALSE]
  colnames(tmp.data ) =  tmp.snpname1
  result.current <- permute.label(label= Label.sub, n = M , # Required number of permutations for a permutation-based p-value 
                                  data= tmp.data , seed =111 )
  p.list.tmp[[tmp.listname]] <- result.current
  p.list.tmp
  
}

stopCluster(cl)


## 9.3 Summarize permuation results -------- 
# library(stringr)
length(unlist(single.p.list))  
tmp <- names(unlist(single.p.list))
tmp.df <- Reduce(rbind, strsplit(tmp , split=c(" ")))
tmp.df2 <- Reduce(rbind, strsplit(tmp.df[,4] , split= "[.]")) # chr, snp rs

snpname.vec <- c(as.character(tmp.df2[,2]))
chr.vec <- as.numeric(as.character(tmp.df[,2]))
B.vec <-  as.numeric(as.character(tmp.df2[,1]))
pvalue.vec <-unlist(single.p.list)

summary.individual.p.df <- as.data.frame(cbind(snpname.vec, chr.vec, B.vec , pvalue.vec ))
summary.individual.p.df$snpname.vec <- as.character(summary.individual.p.df$snpname.vec)
summary.individual.p.df$chr.vec <- as.numeric(as.character(summary.individual.p.df$chr.vec))
summary.individual.p.df$B.vec  <- as.numeric(as.character(summary.individual.p.df$B.vec ))
summary.individual.p.df$pvalue.vec  <- as.numeric(as.character(summary.individual.p.df$pvalue.vec ))
head(summary.individual.p.df)


# Final summary 
final.sum <- 
  dplyr::right_join(x =iffcond.df , y =summary.individual.p.df ,
                               by = c("snpname" = "snpname.vec", 
                                      "chr"="chr.vec", "B" ="B.vec"))

# Add indicator for selected SNPs
final.sum$select <- 0
final.sum$select5[(final.sum$snpname %in% as.character(select.weak.iff$snpname)) | (final.sum$sw=="strong")]  <- 1
sum(final.sum$select)   

# Significant SNPs
# Permutation p value threshold
p.perm.thred = 0.05/5679 # For example, in our study 
# Determine significance for each SNP
final.select <- final.sum[final.sum$select == 1 , ]
final.select$pvalue.vec[final.select$sw =="weak"]
length(which(final.select$pvalue.vec < p.perm.thred ))   
length(which(final.select$pvalue.vec[final.select$sw =="weak"] < p.perm.thred ))  
length(which(final.select$pvalue.vec[final.select$sw =="strong"] < p.perm.thred ))  

