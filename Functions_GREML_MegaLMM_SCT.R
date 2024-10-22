
Solver_GREML = function(Y,Z){
  G = tcrossprod(Z)
  G = G/mean(diag(G))
  diag(G) = diag(G)+0.01
  rownames(G) = colnames(G) = paste('GE',1:nrow(Y),sep='_')
  k = ncol(Y)
  Udta = data.frame(Y,GE=factor(paste('GE',1:nrow(Y),sep='_')))
  asreml.options(gammaPar = FALSE, ai.sing = TRUE)
  a0 = suppressWarnings(
    asreml(fixed = as.formula(paste('cbind(',paste(colnames(Y),collapse = ','),')~trait')),
           random = as.formula(~us(trait):vm(GE,G)),
           residual = as.formula(~units:diag(trait)),
           data=Udta,workspace='40Gb'))
  vc0 = summary(a0)$varcomp$component[ 1: ((k*(k-1))/2+k) ] 
  names(vc0) = rownames(summary(a0)$varcomp)[ 1: ((k*(k-1))/2+k) ] 
  Vg = matrix(NA,k,k,byrow = T)
  Vg[upper.tri(Vg,T)] = vc0
  Vg[lower.tri(Vg)] = Vg[upper.tri(Vg)]
  Ve = rev(rev(summary(a0)$varcomp$component)[1:k])
  gebv =matrix(fitted.asreml(a0),ncol=k,byrow=T)
  return(list(gebv=gebv,Vg=Vg,Ve=Ve))}

Solve_SCT = function(Y,Z){
  EY = bWGR::EigenBDCSVD(Y)
  Q = Y %*% EY$V
  A = bWGR::FUVBETA(Q,Z)
  b = A %*% t(EY$V)
  h = Z %*% b
  return(list(hat=h,b=b))
}

Solve_MegaLMM = function(Y,X){
  # turn into inputs for MegaLMM
  df = data.frame(ID = 1:nrow(Y))
  K = tcrossprod(scale(X,scale=F))/sum(colMeans(X)/2)
  rownames(K) = colnames(K) = rownames(Y) = df$ID
  # basic MegaLMM code
  run_parameters = MegaLMM_control(
    h2_divisions = 20, burn = 0, K = 15)
  MegaLMM_state = setup_model_MegaLMM(
    Y = Y, formula = ~ 1 + (1|ID),
    data = df, relmat = list(ID = K), 
    run_parameters=run_parameters,
    run_ID = 'MegaLMM_basic')
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1 = list(shape = 20,rate=1),
    delta_2 = list(shape = 3,rate=1),
    delta_iterations_factor = 100)
  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 5),   
    tot_F_var = list(V = 18/20, nu = 20),  
    h2_priors_resids_fun = function(h2s,n)  1,  
    h2_priors_factors_fun = function(h2s,n) 1, 
    Lambda_prior = Lambda_prior)
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(Y)+1,verbose=F)
  mapID = which(maps$map_results$total_kept_NAs/length(Y)<0.01)[1]
  MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[mapID]])
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state$run_parameters$burn = run_parameters$burn
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)
  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec','B1')
  MegaLMM_state$Posterior$posteriorMean_params = 'Eta_mean'
  MegaLMM_state$Posterior$posteriorFunctions = list(  U = 'U_F %*% Lambda + U_R + X1 %*% B1',
                                                      h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'  )
  MegaLMM_state = clear_Posterior(MegaLMM_state) 
  n_iter = 100
  MegaLMM_state$run_parameters$which_sampler$Y = 2
  for(i in 1:5) {
    print(sprintf('Burnin run %d',i))
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) 
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)
    traceplot_array(MegaLMM_state$Posterior$Lambda,name = file.path(MegaLMM_state$run_ID,'Lambda.pdf'))
    print(sprintf('Completed %d burnin samples', MegaLMM_state$current_state$nrun))  }
  MegaLMM_state = clear_Posterior(MegaLMM_state)
  n_iter = 250
  for(i in 1:4) {
    print(sprintf('Sampling run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)
    print(MegaLMM_state)  }
  U_samples = load_posterior_param(MegaLMM_state,'U')
  U_hat = get_posterior_mean(U_samples)
  Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
  return(list(gebv=U_hat,hat=Eta_mean))}
