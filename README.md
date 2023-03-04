# mLDA Pipeline for Weak Signal Detection
This repository contains the analysis pipeline of detecting marginally strong and marginally weak signals for binary outcome classification using single nucleotide polymorphism (SNP) genotype data. The developed analytical pipeline is open-source, and can be used for analyzing SNP data in R softwares with step-by-step instructions.


# Overview of Workflow and Pipeline
0) Housekeeping
1) Load outcome labels and demographic data of study population
2) Filter SNP genotype data using minor allele frequency threshold
3) Select marginally strong SNPs using the filtered genotype data
4) Search for marginally strong SNP correlated SNPs based on LD
5) Construct SNP clusters based on LD matrices
6) Prepare data which contains all SNPs in the detected SNP clusters 
7) Detect marginally weak signals based on cluster-adjusted effect sizes
8) Binary outcome classification based on selected SNPs
9) Permutation p value for each individual SNP


# Shorthands: 
Single nucleotide polymorphism (SNP)

Chromosome (chr)

Minor allele frequency (MAF)

Linkage disequilibrium (LD)
