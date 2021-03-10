###
### Files from https://github.com/darenwang/vectordp
### 6-162: dp_funcitons.R
### 181-343: localscreen.R
### 356-468: main.R

####
####
#### dp_functions.R
####
####
#compute res given interval [s,e]

residualvar=function(X, Y,s,e,p,lambda){
  options(warn=-1)
  estimate=matrix(0,nrow=p,ncol=p)
  for ( m in 1:p){
    out=glmnet::glmnet(x=t(X[,  s  : e ]), y=Y[m,  s  :e], family=c("gaussian"),
                       alpha = 1, lambda=lambda/sqrt(e-s) )#,intercept=F)
    estimate[m, ] = as.vector(out$beta)
  }
  y.estimate = estimate%*%X[, s:e ]
  return(norm(y.estimate -Y[, s:e ], type="F" )^2)
}
###end of function



# pull the change point given partitions
dp.pull.change= function(N,  partition){
  cc =N
  estimate.change=cc
  while(cc[1]!=1){
    estimate.change=c(partition[cc],estimate.change)
    cc=partition[cc]
  }
  return(estimate.change[-length(estimate.change)]-1)

}
###end of function

#dp matrix  function
dp.matrix.function=function(X.train, Y.train, p, lambda.lasso , delta1, delta2){
  N=ncol(X.train)
  dp.matrix=matrix(Inf, N, N)
  for (s in 1:N){
    for (e in 1:N){
      if (e - s > delta1 && e - s <= delta2){
        #print(c(i,j))
        dp.matrix[s,e]=residualvar(X=X.train, Y=Y.train,s ,e ,p,lambda.lasso)
      }
    }
  }
  return(dp.matrix)
}

#end

#compute res plus gam
res.and.gam=function(r, l , gam, Best.value, dp.matrix){
  result = ifelse(l == 1, - gam, Best.value[l]) + gam + dp.matrix[l,r]
  return(result)
}

###end of function


dp.gam.function=function(dp.matrix, gam, delta1){
  N=ncol(dp.matrix)
  Best.value=rep(Inf, N)
  partition= rep(1,N )

  for( r in 2:N){
    b=sapply(1:r,function(l)
      res.and.gam(r,l, gam,Best.value,dp.matrix  )
    )
    if (min(b) <Best.value[r]){Best.value[r]=min(b)
    partition[r]=which.min(b)}
  }

  dp.estimate=dp.pull.change(N, partition ) #pull out change points
  return(c( dp.estimate ,N))
}

dp.gam.list=function(X.train, Y.train, p, lambda.lasso, gam.list, delta1, delta2) {
  dp.matrix = dp.matrix.function (X.train, Y.train, p, lambda.lasso , delta1, delta2)
  dp.list = vector('list',length= length(gam.list))
  for ( i in 1: length(gam.list)){
    dp.list[[i]] = dp.gam.function(dp.matrix ,gam.list[i], delta1)
  }
  return(dp.list)
}

compute.test.errors=function(X.train, Y.train,  X.test,Y.test,lambda.lasso,estimate.changes){
  p <- nrow(X.train)
  res=0
  for ( l in 2:length(estimate.changes)){

    s=estimate.changes[l-1]+1
    e=estimate.changes[l]
    estimate=matrix(0,nrow=p,ncol=p)
    for(m in 1:p) {
      out=glmnet::glmnet(x=t(X.train[,  s  : e ]), y=Y.train[m,  s  :e], family=c("gaussian"),
                 alpha = 1, lambda=lambda.lasso/sqrt(e-s )) #,intercept=F)
      estimate[m, ]=as.vector(out$beta )
    }
    res=res+ norm( Y.test[,s  :e] -estimate%*%X.test[, s:e ], type="F" )^2

  }
  return(res)
}

dp.main=function(X.train, Y.train, X.test, Y.test, p, lambda.lasso.list, gam.list, delta1, delta2){
  rec.temp=rep(0, nrow=length(lambda.lasso.list) )
  dp.lambda.list=vector('list', length=length(lambda.lasso.list))
  for (ll in 1:length(lambda.lasso.list)){
    #print(ll)
    dp.list= dp.gam.list (X.train, Y.train, p, lambda.lasso.list[ll] ,gam.list, delta1, delta2)
    res.vec.temp= sapply(1:length(gam.list), function(j)
      compute.test.errors (X.train, Y.train,  X.test,Y.test,lambda.lasso.list[[ll]],dp.list[[j]]) )
    rec.temp[ll]=min(res.vec.temp)
    dp.lambda.list[[ll]]= dp.list[[which.min(res.vec.temp)]]
  }


  return(dp.lambda.list[[which.min( rec.temp)]])

}

#end of dp main


#use in the simulations
dp.main.function=function(data, p, lambda.lasso.list, gam.list, delta1, delta2){
  data.temp=data
  if (ncol(data)%%2 ==0){
    data.temp=data[,2:ncol(data)]
  }

  X=data.temp[,1:(ncol(data.temp) - 1)]
  Y=data.temp[,2:ncol(data.temp)]
  X.train = X[, seq(1, ncol(X), 2)]
  X.test  = X[, seq(1, ncol(X) - 1, 2) + 1]
  Y.train = Y[, seq(1, ncol(Y), 2)]
  Y.test  = Y[, seq(1, ncol(Y) - 1, 2) + 1]
  return(2*dp.main(X.train, Y.train,X.test,Y.test, p, lambda.lasso.list ,gam.list, delta1, delta2))

}












####
####
#### localscreen.R
####
####
#one change point


convert.design.matrix.one.change=function(X,t){
  p <- nrow(X)
  ee=ncol(X)
  xx1=t(X)
  xx1[ (t+1):ee,]=0

  xx2=t(X)
  xx2[1:t,]=0

  xx=cbind(xx1/sqrt(ee),xx2/sqrt(ee))
  index=c()
  for ( pp in 1:p){
    index=c(index,pp,pp+p)
  }
  xxx=xx[,index]
  return(xxx)

}

train.res.group.lasso =function(t, y, X, lambda.group){
  p <- nrow(X)
  group=rep(1:p,each=2)
  out=gglasso::gglasso(x=convert.design.matrix.one.change(X,t),y=y,group=group, loss="ls",
                       lambda=lambda.group/ncol(X) ,intercept = FALSE,eps = 0.0001)
  alpha=as.vector(out$beta)
  #test
  #beta=c(alpha[1:p ]/sqrt(t),alpha[(p+1):(2*p)]/sqrt(ncol(X)-t))
  #beta
  res=sum(( y - convert.design.matrix.one.change(X,t) %*%alpha  )^2)
  return(res)
}

test.res.group.lasso =function(t,y, X, y.test, X.test   ,lambda.group){
  p <- nrow(X)
  group=rep(1:p,each=2)
  out=gglasso::gglasso(x=convert.design.matrix.one.change(X,t),y=y,group=group, loss="ls",
                       lambda=lambda.group/ ncol(X) ,intercept = FALSE,eps = 0.001)
  #test
  alpha=as.vector(out$beta)
  beta=c(alpha[1:p ]/sqrt(t),alpha[(p+1):(2*p)]/sqrt(ncol(X)-t))
  beta
  res=sum(( y.test - convert.design.matrix.one.change(X.test,t) %*%alpha  )^2)
  return(res)
}



res.at.t=function(s,e,t,X.train,Y.train,lambda.group ){
  p <- nrow(X.train)
  return(sum(sapply( 1:p, function(m)
    train.res.group.lasso(t , y=Y.train[m,   s :e], X = X.train[, s: e ]  ,lambda.group ) )))
}


find.one.change.grouplasso=function(s,e,X.train,Y.train,delta.local   ,lambda.group ){
  p <- nrow(X.train)

  estimate= (s+e)/2

  if( e-s > 2*delta.local){
    can.vec=c((s+delta.local): (e-delta.local))
    can.vec=can.vec[which(can.vec%%  2==0)]
    res.seq=sapply(can.vec,function(t)
      res.at.t(s,e,t-s+2  ,X.train,Y.train,lambda.group))
    #plot(can.vec,res.seq)
    estimate= can.vec[which.min(res.seq)]
    #lv= min(res.seq)
  }

  return(  estimate )

}


group.lasso.test.error=function(s,e,X.train, Y.train, X.test, Y.test, delta.local, lambda.group.list ){
  p <- nrow(X.train)
  estimate.temp=rep(0,length(lambda.group.list) )
  test.errors.temp=rep(0,length(lambda.group.list) )
  for ( ll in 1:length(lambda.group.list)){

    estimate.temp[ll]=find.one.change.grouplasso(s,e,X.train,Y.train,delta.local,lambda.group.list[ll] )
    test.errors.temp[ll]=sum(sapply( 1:p, function(m)
      test.res.group.lasso( estimate.temp[ll]-s+2 ,y=Y.train[m,s:e],X= X.train[, s: e  ],
                            y.test=Y.test[m, s:e], X.test= X.test[, s: e  ],lambda.group.list[ll] ) ))
    #print(estimate.temp[ll])
  }

  return(lambda.group.list[which.min(test.errors.temp)])
}

local.group=function(X, Y,X.train, Y.train, X.test, Y.test,delta.local, lambda.group.list,estimated.changes ){
  p <- nrow(X.train)
  result.changes= estimated.changes

  KK=length(estimated.changes)-2
  if (KK>0){
    result.changes[1]=1
    for ( kk in 1:KK){
      #print(kk)
      lambda.temp=group.lasso.test.error (s=result.changes[kk],e=result.changes[kk+2],X.train, Y.train, X.test, Y.test, delta.local, lambda.group.list )
      temp.estimate=
        find.one.change.grouplasso(s=2*result.changes[kk],e=2*result.changes[kk+2],X ,Y ,delta.local   ,lambda.group= lambda.temp )
      result.changes[kk+1]=round(temp.estimate/2)
      #print(temp.estimate)
    }
  }
  result.changes[1]=0
  return(result.changes)

}


local.main.function= function (data,delta.local, lambda.group.list,dp.estimate){
  data.temp=data
  if (ncol(data)%%2 ==0){
    data.temp=data[,2:ncol(data)]
  }

  temp.estimate=dp.estimate/2

  X=data.temp[,1:(ncol(data.temp)-1)]
  Y=data.temp[,2:ncol(data.temp)]
  X.train= X[,seq(1,ncol( X),2)]
  X.test= X[,seq(1,ncol( X)-1,2)+1]
  Y.train= Y[,seq(1,ncol( Y),2)]
  Y.test= Y[,seq(1,ncol( Y)-1,2)+1]
  return(2* local.group (X, Y,X.train, Y.train, X.test, Y.test,delta.local, lambda.group.list,temp.estimate))
}













####
####
#### main.R
####
####

#' @export
var_pgl <- function(x,
                    gamma        = seq(1, 50, 2),
                    lambda       = 0.1 * sqrt(log(p)),
                    lambda_group = 0.3 * sqrt(log(p)),
                    minsl        = 5,
                    maxsl        = round(n / 3)) {
  p <- ncol(x)
  n <- nrow(x)
  # n <- nrow(x)
  x <- t(x)

  # lambda.group.list <- c(0.5, 1, 1.5)
  # minimal sample size for lasso and group lasso estiamtors; needed for any lasso estimators
  # delta1 <- minsl
  # delta.local <- minsl

  # start.time <- Sys.time()
  dp.estimate <- dp.main.function(x, p, lambda, gamma, minsl, maxsl)
  # end.time <- Sys.time();time.taken <- end.time - start.time;print(time.taken)
  lr.estimate <- local.main.function(x, minsl, lambda_group, dp.estimate)
  # 1st and last item among estimates are 1 and n.
  lr.estimate[-c(1, length(lr.estimate))]
}

# p = 10, n = 1000
# > microbenchmark(var_pgl(test_x, maxsl = 1000), var_pgl(test_x, maxsl = 100), times = 1)
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# var_pgl(test_x, maxsl = 1000) 953.5667 953.5667 953.5667 953.5667 953.5667 953.5667     1
# var_pgl(test_x, maxsl = 100) 364.1470 364.1470 364.1470 364.1470 364.1470 364.1470     1
# > 363 / 60
# [1] 6.05
# > 363 / 60 * 100
# [1] 605
# > 363 * 100 / 60 / 60
# [1] 10.08333


# p = 10, n = 100
# > microbenchmark(var_pgl(test_x, maxsl = 100), times = 2)
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# var_pgl(test_x, maxsl = 100) 9.821605 9.821605 10.00587 10.00587 10.19015 10.19015     2
# > 10 * 100 / 60 / 60
# [1] 0.2777778
