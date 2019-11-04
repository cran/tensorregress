#' Simulation of tensor regression models
#'
#' Generate response tensors with multiple covariates under different simulation models, specifically for tensors with 3 modes
#' @param seed    seed for generating data
#' @param whole_shape  a vector containing dimension of the tensor
#' @param core_shape   a vector containing Tucker rank of the coefficient tensor
#' @param p            a vector containing numbers of covariates on each mode, see "details"
#' @param dist         distribution of response tensor, see "details"
#' @param dup          number of simulated tensors from the same linear predictor
#' @param signal       a scalar controlling the max norm of the linear predictor
#' @param block        a vector containing boolean variables, see "details"
#' @return     a list containing the following:
#'
#' \code{tsr} {a list of simulated tensors, with the number of replicates specified by \code{dup}}
#'
#' \code{X_covar1}  {a matrix, covariate on first mode}
#'
#' \code{X_covar2}  {a matrix, covariate on second mode}
#'
#' \code{X_covar3}  {a matrix, covariate on third mode}
#'
#' \code{W} {a list of orthogonal coefficient matrices - one for each mode, with the number of columns given by \code{core_shape}}
#'
#' \code{G}  {an array, core tensor with size specified by \code{core_shape}}
#'
#' \code{C_ts}  {an array, coefficient tensor, Tucker product of \code{G},\code{A},\code{B},\code{C}}
#'
#' \code{U} {an array, linear predictor,i.e. Tucker product of \code{C_ts},\code{X_covar1},\code{X_covar2},\code{X_covar3}}
#'
#' @details    By default non-positive entry in \code{p} indicates no covariate on the corresponding mode of the tensor.
#'
#'             \code{dist} specifies three distributions of response tensor: binary, poisson and normal distribution.
#'
#'             \code{block} specifies whether coefficient factor matrix is a membership matrix, set to \code{TRUE} when utilizing the stochastic block model
#'
#'
#' @export
#' @examples 
#' seed = 34
#' dist = 'binary'
#' data=sim_data(seed, whole_shape = c(20,20,20), core_shape=c(3,3,3),
#' p=c(5,5,5),dist=dist, dup=5, signal=4)


###----  functions for simulation
#####---- This is the function used for generating data through different distribution
#         of core tensor in  semi-supervised setting
## p is the dimension of the covaraite. p = 0 represents the case without covaraites
sim_data = function(seed, whole_shape = c(20,20,20), core_shape = c(3,3,3),p=c(3,3,0),dist, dup, signal,block=rep(FALSE,3)){

  d1 = whole_shape[1] ; d2 = whole_shape[2] ; d3 = whole_shape[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  p1 = p[1]; p2 = p[2]; p3 = p[3];
  
  #### warning for r should larger than 0
  if(r1<=0 | r2 <= 0|r3<= 0){
    warning("the rank of coefficient tensor should be larger than 0",immediate. = T)
    return()
  }
  
  if(p1 > d1|p2 > d2|p3 > d3){
    warning("the number of covariates at each mode should be no larger than the dimension of the tensor",immediate. = T)
  }
  
  #### warning for p should larger than r
  if((p1<r1 & p1>0)|(p2<r2 & p2>0)|(p3<r3 & p3>0)){
    warning("the rank of coefficient tensor should be no larger than the number of covariates",immediate. = T)
  }
  
  
  
  ####-------- generate data
  set.seed(seed)  # 24 # 37  #  347
  X_covar1 = X_covar2 = X_covar3 = NULL

  if(p1<=0){
    X_covar1=diag(1,d1)
    p1=d1
  }else{
    X_covar1 = matrix(rnorm(d1*p1,mean = 0, sd = 1/sqrt(d1)),d1,p1)
  }

  if(p2<=0){
    X_covar2=diag(1,d2)
    p2=d2
  }else{
    X_covar2 = matrix(rnorm(d2*p2,mean = 0, sd =1/sqrt(d2)),d2,p2)
  }


  if(p3<=0){
    X_covar3=diag(1,d3)
    p3=d3
  }else{
    X_covar3 = matrix(rnorm(d3*p3,mean = 0, sd = 1/sqrt(d3)),d3,p3)
  }
  
  #### if use block, p and r should larger than 1 of smaller than 0
  if(block[1] == T & (p1 ==1|r1 == 1)){
    warning("number of groups should be larger than 1",immediate. = T)
  }
  
  if(block[2] == T & (p2 ==1|r2 == 1)){
    warning("number of groups should be larger than 1",immediate. = T)
  }
  
  if(block[3] == T & (p3 ==1|r3 == 1)){
    warning("number of groups should be larger than 1",immediate. = T)
  }
  
  
  if (block[1]==TRUE){
    b1=sort(sample(1:r1,p1,replace=TRUE))
    if(length(unique(b1) )== 1){
      warning("mode-1 degenerates to 1",immediate. = T)
    }
    W1=model.matrix(~-1+as.factor(b1))
    r1=dim(W1)[2]
  }else W1 =as.matrix(randortho(p1)[,1:r1])

  A = X_covar1%*%W1 ## factor matrix

  if (block[2]==TRUE){
    b2=sort(sample(1:r2,p2,replace=TRUE))
    if(length(unique(b2)) == 1){
      warning("mode-2 degenerates to 1",immediate. = T)
    }
    W2=model.matrix(~-1+as.factor(b2))
    r2=dim(W2)[2]
  }else W2 = as.matrix(randortho(p2)[,1:r2])
  B = X_covar2%*%W2 ## factor matrix

  if (block[3]==TRUE){
    b3=sort(sample(1:r3,p3,replace=TRUE))
    if(length(unique(b3)) == 1){
      warning("mode-3 degenerates to 1",immediate. = T)
    }
    W3=model.matrix(~-1+as.factor(b3))
    r3=dim(W3)[2]
  }else W3 = as.matrix(randortho(p3)[,1:r3])
  C= X_covar3%*%W3 ## factor matrix
  

  ### G: core tensor
  G = as.tensor(array(runif(r1*r2*r3,min=-1,max=1),dim = c(r1,r2,r3)))


  ### U: linear predictor
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  G=G/max(abs(U))*signal ## rescale subject to entrywise constraint
  U=U/max(abs(U))*signal

  C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data ## coefficient

  ### tsr:binary tensor
  if(dist=="binary"){
    tsr = lapply(seq(dup), function(x) array(rbinom(d1*d2*d3,1,prob = as.vector( 1/(1 + exp(-U)))),dim = c(d1,d2,d3)))}
  else if (dist=="normal"){
    tsr = lapply(seq(dup), function(x) array(rnorm(d1*d2*d3,U),dim = c(d1,d2,d3)))}#sd = 1
  else if (dist=="poisson"){
    tsr = lapply(seq(dup), function(x) array(rpois(d1*d2*d3,exp(U)),dim = c(d1,d2,d3)))}


  return(list(tsr = tsr,X_covar1 = X_covar1, X_covar2 = X_covar2,X_covar3 = X_covar3,
              W = list(W1 = W1,W2 = W2,W3 = W3), G=G@data, U=U,C_ts=C_ts))
}









