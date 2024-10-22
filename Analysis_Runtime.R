
# Load packages
require(bWGR)
require(MegaLMM)
require(asreml)

# Simulate genomic data
Z = SimZ(ind = 500, # number of individuals
         chr = 10,   # number of chromosomes
         snp = 150, # number of markers per chromosomes
         F2 = TRUE)   # simulate a F2 population

# Load MegaLMM function
source('Functions_GREML_MegaLMM_SCT.R')

# Generate Qa parameterization
V = EigenBDCSVD(EigenCNT(Z))$V;
Q = EigenCNT(Z) %*% V;

######## SIMULATE A SCENARIO #############

# Data sparseness
seed = 123
set.seed(seed)
MissPer = 0

# Simulate GxE correlation
num_of_environments = 100
GC = SimGC(k = num_of_environments)
image(GC)

# Set heritability value
h2 = 0.2

# Simulate phenotypes
S = SimY(Z = Z, 
         GC = GC, 
         h2 = h2,
         seed = seed)

# Extract matrix of phenotypes and TBV
Y = S$Y
TBV = S$tbv

# UV
uv_time = system.time(fx<-apply(Y,2,function(y) MRR3(matrix(y),Q)$hat))[3]
uv_acc = sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(uv_time,'/',round(mean(uv_acc),3),'\n')

# SCT (simplified canonical transformation)
sct_time = system.time(fx<-Solve_SCT(Y,Q)$hat)[3]
sct_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(sct_time,'/',round(mean(sct_acc),3),'\n')

# GREML
greml_time = system.time(fx<-Solver_GREML(Y,Q)$gebv)[3]
greml_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(greml_time,'/',round(mean(greml_acc),3),'\n')

# HCS-PEGS
hcs_time = system.time(fx<-MRR3(Y,Q,HCS=T,cores=-1,tol=1e-09)$hat)[3]
hcs_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(hcs_time,'/',round(mean(hcs_acc),3),'\n')

# XFA-PEGS
xfa_time = system.time(fx<-MRR3(Y,Q,XFA=T,TH=T,NumXFA=3,cores=-1,tol=1e-09)$hat)[3]
if(anyNA(fx)){cat('fallback ');rm(fx);gc(); xfa_time = system.time(fx<-MRR3F(Y,Q,XFA=T,TH=T,NumXFA=3)$hat)[3];}
xfa_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(xfa_time,'/',round(mean(xfa_acc),3),'\n')

# PEGS (unstructured covariance)
mv_time = system.time(fx<-MRR3(Y,Q,TH=T,cores=-1,tol=1e-09)$hat)[3]
mv_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(mv_time,'/',round(mean(mv_acc),3),'\n')

# MegaSEM (Full-rank model)
megasem_time = system.time(fx<-ZSEMF(Y,Q)$hat)[3]
megasem_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(megasem_time,'/',round(mean(megasem_acc),3),'\n')

# MegaLMM
megalmm_time = system.time(fx <- try(Solve_MegaLMM(Y,Z)$gebv,silent = T))[3]
megalmm_acc =  sapply(1:num_of_environments, function(i) cor(TBV[,i],fx[,i]) )
cat(megalmm_time,'/',round(mean(megasem_lmm),3),'\n')

# to lower runtime of bWGR functions, replace MRR3 by MRR3F
