
# Load packages
require(bWGR)
require(MegaLMM)

# Load data
load('SoyNAM_data.RData')

# Generate Qa parameterization
V = EigenBDCSVD(EigenCNT(X))$V;
Q = EigenCNT(X) %*% V;

# Load MegaLMM function
source('GREML_MegaLMM_SCT_functions.R')

######## SIMULATE A SCENARIO ##############

# Data sparseness
seed = 123
set.seed(seed)
MissPer = 0.75

# Simulate GxE correlations
GC = SimGC(100)
image(GC)
hist(GC)

# Set heritability value
h2 = 0.4

# Simulate phenotype with unstructured 
S = SimY(X,
         k=100, 
         GC=GC, # input value to fix GxE, e.g. GC=0.5
         h2=h2,
         PercMiss=MissPer,
         seed=seed,
         BlkMiss=TRUE)

# Extract matrix of phenotypes and TBV
Y = S$Y
TBV = S$tbv

# UV
uv_time = system.time(fx<-apply(Y,2,function(y) MRR3(matrix(y),Q)$hat))[3]
uv_acc = colMeans(sapply(1:ncol(TBV), function(i) by(data.frame(fx[,i],TBV[,i]),pop,function(x) cor(x)[1,2] ))) 
cat(uv_time,'/',round(mean(uv_acc),3),'\n')

# XFA-PEGS
xfa_time = system.time(fx<-MRR3(Y,Q,XFA=T,TH=T,NumXFA=3,cores=-1,tol=1e-09))[3]
if(anyNA(fx$hat)){cat('fallback ');rm(fx);gc(); xfa_time = system.time(fx<-MRR3F(Y,Q,XFA=T,TH=T,NumXFA=3))[3];}
xfa_acc = colMeans(sapply(1:ncol(TBV), function(i) by(data.frame(fx$hat[,i],TBV[,i]),pop,function(x) cor(x)[1,2] ))) 
cat(xfa_time,'/',round(mean(xfa_acc),3),'\n')

# HCS-PEGS
cat('HCS-PEGS ')
hcs_time = system.time(fx<-MRR3(Y,Q,HCS=T,cores=-1,tol=1e-09))[3]
hcs_acc = colMeans(sapply(1:ncol(TBV), function(i) by(data.frame(fx$hat[,i],TBV[,i]),pop,function(x) cor(x)[1,2] ))) 
cat(hcs_time,'/',round(mean(hcs_acc),3),'\n')

# # PEGS
mv_time = system.time(fx<-MRR3(Y,Q,TH=T,cores=-1,tol=1e-09))[3]
mv_acc = colMeans(sapply(1:ncol(TBV), function(i) by(data.frame(fx$hat[,i],TBV[,i]),pop,function(x) cor(x)[1,2] ))) 
cat(mv_time,'/',round(mean(mv_acc),3),'\n')

# MegaSEM
megasem_time = system.time(fx<-ZSEMF(Y,Q))[3]
megasem_acc = colMeans(sapply(1:ncol(TBV), function(i) by(data.frame(fx$hat[,i],TBV[,i]),pop,function(x) cor(x)[1,2] ))) 
cat(megasem_time,'/',round(mean(megasem_acc),3),'\n')

# MegaLMM
megalmm_time = system.time(fx <- try(Solve_MegaLMM(Y,X),silent = T))[3]
megalmm_acc = colMeans(sapply(1:ncol(TBV), function(i) by(data.frame(fx$gebv[,i],TBV[,i]),pop,function(x) cor(x)[1,2] ))) 
cat(megalmm_time,'/',round(mean(megasem_lmm),3),'\n')
