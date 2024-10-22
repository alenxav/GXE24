
# Load packages
require(bWGR)
require(MegaLMM)

# Load data
load("Data_G2F.RData")

# Load MegaLMM, SCT, GREML functions
source('Functions_GREML_MegaLMM_SCT.R')

# Generate Qa parameterization
K = EigenGRM(Z)
system.time(Q<-K2X(K))[3]
Q = Q[,which(apply(Q,2,var)>1e-6)] # drop eigen pairs with low variance
rownames(Q) = rownames(Z)

# Individuals in the estimation and prediction set
es = intersect(rownames(Y_ES),rownames(Q))
ps = intersect(rownames(Y_PS),rownames(Q))

# Validation function (summarize results)
validation = function(hat,avg=TRUE){
  colnames(hat) = colnames(Y_ES)
  st1 = substr(colnames(hat),1,4)
  st2 = substr(colnames(Y_PS),1,4)
  prd = t(apply(hat,1,function(x) tapply(x,st1,mean))[st2,])
  x = data.frame(
    ByState = diag(cor(prd,Y_PS,use='p')),
    AllxAll =  colMeans(cor(hat,Y_PS,use='p')),
    Average =  cor(rowMeans(hat),Y_PS,use='p')[1,])
  xx = apply(x,2,function(xx){
    xx = round(c(Mu=mean(xx),Std=sd(xx)),2)
    paste0(xx[1],' (',xx[2],')')}) 
  if(avg) return(xx) else return(x)
}

######## Train and validate models ##########

# UVW (Univariate computed within-environment)
beta_uvw = FUVBETA(Y_ES[es,],Q[es,])
uvw_pred = Q[ps,] %*% beta_uvw
validation(uvw_pred)

# MegaSEM (Full-rank model)
beta_sem =  ZSEMF(Y_ES[es,],Q[es,])
pred_sem = Q[ps,] %*% beta_sem$b
validation(pred_sem)

# UVA (Univariate computed across environments)
tall = reshape2::melt(Y_ES,na.rm = T)
blup = mixed(value,~Var1,~Var2,tall)$Coefficients$Var1
genotyped_individuals = intersect(names(blup),rownames(Q))
beta_uva = MRR3F(Y = as.matrix(blup[genotyped_individuals]),
                 X = Q[genotyped_individuals,])$b
uva_pred = c(Q[ps,] %*% beta_uva) 
uva_result = cor(uva_pred,Y_PS,use='p')
round(c(Mu=mean(uva_result),Std=sd(uva_result)),2)

# HCS (Multivariate heteroskedastic compound symmetry structured covariance)
beta_hcs =  MRR3F(Y_ES[es,],Q[es,],HCS=TRUE,TH=TRUE,cores=-1)
pred_hcs = Q[ps,] %*% beta_hcs$b
validation(pred_hcs)

# MV (Multivariate unstructured covariance)
beta_mv =  MRR3F(Y_ES[es,],Q[es,],cores=-1,TH=TRUE)
pred_mv = Q[ps,] %*% beta_mv$b
validation(pred_mv)

# XFA (Multivariate extended factor analytics structured covariance)
beta_xfa =  MRR3F(Y_ES[es,],Q[es,],XFA=TRUE,TH=TRUE,cores=-1)
pred_xfa = Q[ps,] %*% beta_xfa$b
validation(pred_xfa)
