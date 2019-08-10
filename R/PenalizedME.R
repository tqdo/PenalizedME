# Estimate the optimal coefficients for the Penalized mixture of experts (Do & Gregory, 2019+).
#
# @param x the design matrix which does not contain the intercept
# @param y the response vector (continuous only)
# @param lamb.max the maximum value of lambda in the sequence
# @param len the number of different values of lambdas to be included in the sequence
# This function returns a list which contains estimated coefficients for gating networks and experts
PenalizedME<- function(x,y,lamb.max=2.5,len=20){
  library(Matrix)
  library(MASS)
  library(NominalLogisticBiplot)
  library(gdata)
  library(stats)
  library(glmnet)
  eps = 10^(-5)
  n1 = nrow(x)
  set.seed(1)
  b = cbind(x,y)[sample(n1),]
  x1 = x
  y1 = y
  train.index1 = 1:(round(.7*n1,0))
  valid.index1 = (max(train.index1)+1):(max(train.index1)+round((n1-max(train.index1))/2,0))
  test.index1 = (max(valid.index1)+1):n1

  x1train = x1[train.index1,]
  x1valid = x1[valid.index1,]
  x1test = x1[test.index1,]

  y1train = y1[train.index1]
  y1valid = y1[valid.index1]
  y1test = y1[test.index1]
  sample=c(train.index5,valid.index5)
  unscale.di=b[sample,]
  mean.v=apply(unscale.di,2,mean)
  sd.v=apply(unscale.di,2,sd)
  di = scale(unscale.di)
  y=di[,ncol(di)] # training response
  x=cbind(rep(1,length(y)),di[,1:(ncol(di)-1)]) # training predictors
  dim = ncol(x) # dimensions
  n = nrow(x) # no of observations
  no.exp = 6   # number of experts
  test = (b[-sample,]-matrix(rep(mean.v,nrow(b)-length(sample)),ncol=dim,byrow=T))/matrix(rep(sd.v,nrow(b)-length(sample)),ncol=dim,byrow=T)

  lambda=seq(0.0,lamb.max,length.out = len)
  # Initial parameters
  l=0.1 # a constant
  sd=10 # sigma

  gn.para = matrix(c(rep(0.1,dim*(no.exp-1)),rep(0,dim)),nrow =dim,ncol=no.exp) # Parameters for gating networks
  exp.para = matrix(l,nrow =dim,ncol=no.exp ) # Parameters for experts
  set.seed(1)

  km = kmeans(x[,-1],no.exp)
  cl = km$cluster
  ridgebetas = solve(t(x)%*%x+0.0001*diag(rep(1,dim)))%*%t(x)%*%y
  for(i in 1:no.exp){
    exp.para[,i] = solve(t(x[cl==i,])%*%x[cl==i,]+0.0001*diag(rep(1,dim)))%*%t(x[cl==i,])%*%y[cl==i]
    exp.para[,i][abs(exp.para[,i])<0.0001] = ridgebetas[which(abs(exp.para[,i])<0.0001)]
  }
  new.exp.para = exp.para # Just some logistic stuff
  new.gn.para = gn.para # Just some logistic stuff

  #Calculate the weight (posterior)

  #Gating network
  p.i = exp(as.matrix(x) %*% as.matrix(gn.para))
  p.i = p.i/rowSums(p.i)
  p.i.nr = p.i
  #Expert
  p.exp.i = apply(exp.para,2, function(z) dnorm(y,mean = as.matrix(x) %*% as.matrix(z), sd= sd))

  #Posterior (multiple gating and expert)
  post = p.exp.i*p.i
  post = post/rowSums(post)

  longy = as.vector(post[,-no.exp])
  longp = as.vector(p.i[,-no.exp])
  Xtil = bdiag(replicate(no.exp-1, x, simplify = FALSE))
  W = matrix(0,n*(no.exp-1),n*(no.exp-1))
  for(i in (seq(1,n*(no.exp-1),by=n))){
    for(j in (seq(1,n*(no.exp-1),by=n))){
      if(i==j){W[j:(j+n-1),i:(i+n-1)] = diag(n)*(1/2-1/(2*no.exp))}
      if(!i==j) {W[j:(j+n-1),i:(i+n-1)] = -diag(n)/(2*no.exp)}
    }
  }
  sW =solve(W)
  z = Xtil %*% as.vector(gn.para[,-no.exp]) + sW%*% (longy-longp)
  ev <- eigen(W)
  L <- ev$values
  V <- ev$vectors
  LV= t(sqrt(diag(L))) %*% t(V)
  ztil = LV %*% z
  x2til = LV %*% Xtil
  LVX=x2til


  D = matrix(0,nrow=dim*no.exp*(no.exp-1),ncol=2*dim*no.exp)
  Dindex=1
  Drow=1-2*dim
  repeat{
    j=Dindex+2*dim
    repeat{
      Drow=Drow+2*dim
      D[Drow:(Drow+2*dim-1),Dindex:(Dindex+2*dim-1)]=diag(rep(1,2*dim))
      D[Drow:(Drow+2*dim-1),j:(j+2*dim-1)]=-diag(rep(1,2*dim))

      if(j+2*dim-1==2*dim*no.exp){break}
      j=j+2*dim
    }
    Dindex=Dindex+2*dim
    if(Dindex==2*dim*no.exp+1-2*dim){break}
  }

  Cp1 = matrix(0,nrow=2*dim*no.exp,ncol = dim*no.exp )
  Cp2 = matrix(0,nrow=2*dim*no.exp,ncol = dim*no.exp-dim )
  cpi=1
  cpj=1
  repeat{
    Cp1[cpj:(cpj+dim-1),cpi:(cpi+dim-1)] = diag(rep(1,dim))
    cpi = cpi + dim
    cpj = cpj+2*dim
    if(cpj>2*dim*no.exp-dim){break}
  }

  cpi=1
  cpj=1+dim
  repeat{
    Cp2[cpj:(cpj+dim-1),cpi:(cpi+dim-1)] = diag(rep(1,dim))
    cpi = cpi + dim
    cpj = cpj+2*dim
    if(cpj>2*dim*no.exp-dim){break}
  }
  Cp = cbind(Cp1,Cp2)
  D = D %*% Cp

  ystar = diag(sqrt(unlist(as.data.frame(post)))) %*% rep(as.vector(y),no.exp)/(sqrt(2)*sd)
  xstar = data.matrix(diag(sqrt(unlist(as.data.frame(post)))) %*% bdiag(replicate(no.exp, x, simplify = FALSE)))/(sqrt(2)*sd)

  y2star = rbind(ystar,ztil)
  x2star = bdiag(xstar,x2til)

  Xg = ginv(as.matrix(t(x2star) %*% x2star)) %*% t(x2star) # G-inverse of X
  xX = as.matrix(x2star) %*% as.matrix(Xg)
  ytil = x2star %*% Xg %*% y2star
  Dtil = D %*% Xg

  u.exp.para = ginv(as.matrix(t(Dtil))) %*% (ytil - x2star %*% c(as.vector(exp.para),as.vector(gn.para[,-no.exp]))) # Initial values of u
  new.u.exp.para=u.exp.para


  ##Define functions to minimize

  #Gating network
  fr = function(t){
    total.gn.para = matrix(c(t,rep(0,dim)),nrow =dim,ncol=no.exp)
    test.p.i = exp(as.matrix(x) %*% total.gn.para)
    test.p.i = test.p.i/rowSums(test.p.i)
    ln.test.p.i = log(test.p.i)
    -sum(ln.test.p.i*post)
  }
  max.ite=30
  # Set up for EM
  total.mark = NULL  # mark identical experts (0 means different from others, same number
  #  but not 0 means they are identical)

  no.group=cbind(lambda,rep(no.exp,length(lambda))) # number of groups

  #Logistics stuff
  total.exp.para=NULL
  total.gn.para=NULL
  total.lambda= as.vector(rep(lambda,each=dim))
  train.mse.v = NULL
  test.mse.v = NULL

  #Plot the data
  if(dim==2){plot(x[,-1],y)}

  # Count how many times we do orginal least squares without penalty
  check=0

  up = ginv(as.matrix(x2star))

  #Start
  for(la in 1:length(lambda)){
    it=1
    lamb=lambda[la]
    loss=NULL;loss[it]=-sum(log(rowSums( p.exp.i*p.i)))
    #EM

    repeat{
      fast = NULL
      for(i in seq(1,length(u.exp.para),by=2*dim)){
        fast=rbind(fast,as.matrix(solve(as.matrix(t(t(Dtil)[,i:(i+2*dim-1)])%*%t(Dtil)[,i:(i+2*dim-1)]))%*%t(t(Dtil)[,i:(i+2*dim-1)])))

      }
      fasty=fast%*%ytil
      fastD= NULL
      for(i in seq(1,length(u.exp.para),by=2*dim)){
        fastD=rbind(fastD,as.matrix(fast[i:(i+dim*2-1),])%*%as.matrix(t(Dtil)[,-(i:(i+2*dim-1))]))
      }
      if(lamb>0){
        itnr = 1
        repeat{
          diff=NULL
          itt=0
          repeat{
            for(i in seq(1,length(u.exp.para),by=2*dim)){
              inside=fasty[i:(i+dim*2-1),]-fastD[i:(i+dim*2-1),]%*%u.exp.para[-(i:(i+2*dim-1))]
              no = norm(inside,type="2")
              if(no>lamb){
                new.u.exp.para[i:(i+2*dim-1)] = lamb*inside/no
              } else {
                new.u.exp.para[i:(i+2*dim-1)] = inside
              }
              diff[i:(i+2*dim-1)] = new.u.exp.para[i:(i+2*dim-1)]-u.exp.para[i:(i+2*dim-1)]
              u.exp.para[i:(i+2*dim-1)]=new.u.exp.para[i:(i+2*dim-1)]
            }
            if((max(abs(diff))<0.0000001 ) ){break}
            itt = itt + 1
          }
          # Recoer primal solution
          new.para=matrix(up %*% (ytil-t(Dtil) %*% as.vector(u.exp.para)),nrow =dim)
          new.gn.para = cbind(new.para[,(no.exp+1):(2*no.exp-1)],rep(0,dim))


          p.i.nr1 = exp(as.matrix(x) %*% as.matrix(new.gn.para))
          p.i.nr1 = p.i.nr1/rowSums(p.i.nr1)
          if(max(abs(p.i.nr1-p.i.nr))<0.01){break}
          p.i.nr = p.i.nr1
          longp = as.vector(p.i.nr[,-no.exp])
          z =Xtil %*% as.vector(new.gn.para[,-no.exp]) + sW%*% (longy-longp)
          ztil = LV %*% z

          y2star = rbind(ystar,as.matrix(ztil))
          ytil = xX %*% as.vector(y2star)
          itnr = itnr + 1
          #if(itnr>20){break}
          gn.para = new.gn.para
        }
      }
      if(!lamb==0){new.exp.para = new.para[,1:no.exp]  }

      if(lamb==0){
        new.exp.para=matrix(solve(t(xstar) %*% xstar) %*% t(xstar) %*% ystar,nrow =dim,ncol=no.exp)
        #Optimize for gating network, I just use optim instead of spg
        new.gn.para=matrix(c(optim(unlist(gn.para[,-ncol(gn.para)]),fr,method = "Nelder-Mead")$par,rep(0,dim)),nrow =dim,ncol=no.exp)
      }
      new.sd = sqrt(sum((matrix(y, length(y), no.exp)-as.matrix(x) %*% as.matrix(new.exp.para))^2*post)/n)
      #Gating network
      p.i = exp(as.matrix(x) %*% as.matrix(new.gn.para))
      p.i = p.i/rowSums(p.i)
      print(c(Sys.time()));print(5);print("\n")

      #Expert
      p.exp.i = apply(new.exp.para,2, function(z) dnorm(y,mean = as.matrix(x) %*% as.matrix(z), sd= new.sd))
      print(c(Sys.time()));print(6);print("\n")

      #Posterior (multiple gating and expert)
      post1 = p.exp.i*p.i
      post1= post1/rowSums(post1)
      print(c(Sys.time()));print(7);print("\n")

      if((max(abs(post1-post))<0.01) & (it>15) & (lamb==0)){break}
      if((max(abs(post1-post))<0.01)  & (!lamb==0)){break}
      if(it==100){break}
      post = post1
      ### Recalulate Posterior with new parameters
      longy = as.vector(post[,-no.exp])
      longp = as.vector(p.i[,-no.exp])
      z =Xtil %*% as.vector(new.gn.para[,-no.exp]) + sW%*% (longy-longp)

      ztil = LV %*% z

      ystar = diag(sqrt(unlist(as.data.frame(post)))) %*% rep(as.vector(y),no.exp)/(sqrt(2)*new.sd)
      xstar = data.matrix(diag(sqrt(unlist(as.data.frame(post)))) %*% bdiag(replicate(no.exp, x, simplify = FALSE)))/(sqrt(2)*new.sd)
      y2star = rbind(ystar,as.matrix(ztil))
      x2star = bdiag(xstar,x2til)
      Xg = ginv(as.matrix(t(x2star) %*% x2star)) %*% t(x2star) # G-inverse of X
      xX = x2star %*% Xg
      ytil = xX %*% y2star
      Dtil = D %*% Xg
      up = ginv(as.matrix(x2star))
      ### Parameters update ( the dual parameters update themselves)
      gn.para = new.gn.para
      exp.para = new.exp.para
      sd=new.sd

      it = it +1;cat("iteration=",it)
    }

    gn.para = new.gn.para
    exp.para = new.exp.para
    sd=new.sd

    # Print out the index to keep track
    cat("\n")
    cat(la)
    cat("\n")

    #Prediction in training data
    y1=rowSums((data.matrix(x) %*% new.exp.para)*p.i)
    #Traning MSE
    train.mse = mean((y1*sd.v[dim]-y*sd.v[dim])^2)

    #Set up test data
    x2 = test[,-ncol(test)]
    x2=cbind(rep(1,nrow(as.matrix(x2))),x2)

    # Calculate weight for each expert in testing data
    test.p.i = exp(as.matrix(x2) %*% as.matrix(new.gn.para))
    test.p.i = test.p.i/rowSums(test.p.i)

    # Prediction in test data
    y2=rowSums((data.matrix(x2) %*% new.exp.para)*test.p.i)

    #Test MSE
    test.mse = mean((y2*sd.v[dim]-test[,ncol(test)]*sd.v[dim])^2)

    # Save in the vector
    train.mse.v[la]=train.mse
    test.mse.v[la]=test.mse

    # Mark different and identical experts
    th = 0.01 # threshold to measure difference, smaller than this implies identical
    new.mark = rep(0,no.exp) # vector to store the marks
    for(i in 1:(no.exp-1)){
      for(j in (i+1):no.exp){
        if(max(abs(exp.para[,i]-exp.para[,j]))<=th ){
          if(new.mark[i]==0 & new.mark[j]==0){
            new.mark[i] = min(i,j)
            new.mark[j] = new.mark[i]
          } else{
            new.mark[i] = new.mark[c(i,j)[which.max(c(!new.mark[i]==0,!new.mark[j]==0))]]
            new.mark[j] = new.mark[i]
          }
        }
      }
    }

    total.mark=rbind(total.mark,new.mark) # Save it
    no.group[la,2]=sum(unique(new.mark)>0)+sum(new.mark==0) # Count number of groups
    total.exp.para=rbind(total.exp.para,new.exp.para)  # Save all experts parameters
    total.gn.para=rbind(total.gn.para,new.gn.para) # Save all gating network parameters
  }
  total.mark = cbind(lambda,total.mark) # Attach sequence of lambda to the marking
  index = which.min(test.mse.v)
  chosen.gn = total.gn.para[(dim*index-dim+1):(dim*index),]
  chosen.exp = total.exp.para[(dim*index-dim+1):(dim*index),]
  return(list(chosen.gn,chosen.exp))
}
