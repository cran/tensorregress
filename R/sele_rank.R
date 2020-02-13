#' Rank selection
#'
#' Estimate the Tucker rank of coefficient tensor based on BIC criterion. The choice of BIC
#'  aims to balance between the goodness-of-fit for the data and the degree of freedom in the population model.
#' @param tsr    response tensor with 3 modes
#' @param X_covar1    covariate on first mode
#' @param X_covar2    covariate on second mode
#' @param X_covar3    covariate on third mode
#' @param rank_range  a matrix containing rank candidates on each row
#' @param Nsim        max number of iterations if update does not convergence
#' @param cons        the constraint method, "non" for without constraint, "vanilla" for global scale down at each iteration,
#'
#'                    "penalty" for adding log-barrier penalty to object function.
#' @param lambda      penalty coefficient for "penalty" constraint
#' @param alpha       max norm constraint on linear predictor
#' @param solver      solver for solving object function when using "penalty" constraint, see "details"
#' @param dist        distribution of response tensor, see "details"
#' @return     a list containing the following:
#'
#'                    \code{rank} a vector with selected rank with minimal BIC
#'
#'                    \code{result}  a matrix containing rank candidate and its loglikelihood and BIC on each row

#' @details    For rank selection, recommend using non-constraint version.
#'             
#'            Constraint \code{penalty} adds log-barrier regularizer to
#'            general object function (negative log-likelihood). The main function uses solver in function "optim" to
#'            solve the objective function. The "solver" passes to the argument "method" in function "optim".
#'            
#'             \code{dist} specifies three distributions of response tensor: binary, poisson and normal distributions.
#'
#'
#' @export
#' @examples 
#' seed=24
#' dist="binary"
#' data=sim_data(seed, whole_shape = c(20,20,20),
#' core_shape=c(3,3,3),p=c(5,5,5),dist=dist, dup=5, signal=4)
#' rank_range = rbind(c(3,3,3),c(3,3,2),c(3,2,2),c(2,2,2),c(3,2,3)) 
#' re = sele_rank(data$tsr[[1]],data$X_covar1,data$X_covar2,data$X_covar3,
#'  rank_range = rank_range,Nsim=10,cons = 'non',dist = dist)






sele_rank = function(tsr, X_covar1 = NULL, X_covar2 = NULL, X_covar3 = NULL,rank_range,Nsim=10,cons = 'non', lambda = 0.1, alpha = 1, solver ="CG",dist){
  whole_shape=dim(tsr)                      
  p=rep(0,3)
  if(is.null(X_covar1)) p[1]=whole_shape[1] else p[1]=dim(X_covar1)[2]
  if(is.null(X_covar2)) p[2]=whole_shape[2] else p[2]=dim(X_covar2)[2]
  if(is.null(X_covar3)) p[3]=whole_shape[3] else p[3]=dim(X_covar3)[2]
  
  
  rank_matrix=rank_range
  rank=as.matrix(rank_range)
  
  whole_shape = dim(tsr)
  rank = lapply(1:dim(rank)[1], function(x) rank[x,]) ## turn rank to a list
  upp = lapply(rank, FUN= tensor_regress,tsr = tsr,X_covar1 = X_covar1,X_covar2 = X_covar2,X_covar3 = X_covar3, Nsim = Nsim, cons = cons,lambda = lambda, alpha = alpha, solver = solver,dist=dist)
  
  lglk= unlist(lapply(seq(length(upp)), function(x) tail(upp[[x]]$lglk,1)))
  BIC = unlist(lapply(seq(length(rank)), function(x) (prod(rank[[x]]) + sum((p-rank[[x]])*rank[[x]])) * log(prod(whole_shape))))
  BIC = -2*lglk + BIC
  rank_matrix=cbind(rank_matrix,lglk,BIC)
  
  return(list(rank = rank[[which(BIC == min(BIC))]],result=rank_matrix))
}






