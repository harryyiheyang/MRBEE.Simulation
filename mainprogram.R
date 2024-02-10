for(iter in istar:500){
  
  if((iter-1)%%10==0){
    G <- matrix(rbinom(N * M, 2, 0.3)-0.6, nrow = N, ncol = M) / sqrt(2 * 0.21)
  }
  
  gammaX=MASS::mvrnorm(M,mu=rep(0,4),Sigma=SigmaB);gammaX[,3]=-gammaX[,3];
  gammaX=pnorm(gammaX)*0.22;
  gammau=rnorm(M,0,0.1);
  
  if(CHP==0){
    gammau=gammau*0
  }
  if(CHP==1){
    CHPind=sample(M,M,0.7);
    gammau[CHPind]=0;
    gammaX[-CHPind,]=0
  }
  
  etaX=matrixMultiply(G,gammaX);
  etau=matrixVectorMultiply(G,gammau);etaU=cbind(etau,etau,etau,etau)
  etaX=etaX+0.25*etaU
  
  E=MASS::mvrnorm(N,mu=rep(0,6),SigmaE);E[,4]=-E[,4];
  adj=var(etaX[,1])/hx-var(etaX[,1]);adj=adj*(1+0.25^2-2*0.25*0.5)
  u=matrixVectorMultiply(G,gammau)+E[,1]*sqrt(adj)
  X=matrixMultiply(G,gammaX)+cbind(u,u,-u,u)*0.25+E[,2:5]*sqrt(adj)
  
  gammaa=gammau*0
  if(BUHP==1){gammaa=rnorm(M,0,0.1)}
  if(UUHP==1){gammaa=rnorm(M,0.1,0.2);gammaa[sample(M,M,0.7)]=0}
  if(CHP==1){gammaa=rnorm(M,0.0,0.2);gammaa[CHPind]=0}
  
  y=X%*%theta+matrixVectorMultiply(G,gammaa)+u+E[,6]*sqrt(adj)
  fity=biggwas(y[N2],G[N2,])
  fitx1=biggwas(X[N1,1],G[N1,])
  fitx2=biggwas(X[N1,2],G[N1,])
  fitx3=biggwas(X[N1,3],G[N1,])
  fitx4=biggwas(X[N1,4],G[N1,])
  
  Rxy=cor(cbind(X,y))
  Rxx=Rxy[1:4,1:4]
  rxy=Rxy[1:4,5]*overlap
  Rxy[1:4,5]=Rxy[5,1:4]=rxy
  
  by=fity$est;byse=fity$std
  bx=cbind(fitx1$est,fitx2$est,fitx3$est,fitx4$est)
  bxse=cbind(fitx1$std,fitx2$std,fitx3$std,fitx4$std)
  
  tt=Sys.time()
  fit0=MRBEE.IMRP(by=by,bX=bx,byse=byse,bXse=bxse,Rxy=Rxy,FDR=FDR,var.est=var.est,adjust.method=adjust.method)
  hatb0=fit0$theta
  stdb0=sqrt(diag(fit0$covtheta))
  t0=Sys.time()-tt
  
  mvinput=MendelianRandomization::mr_mvinput(bx=bx,bxse=bxse,by=as.numeric(by),byse=as.numeric(fity$std))
  
  tt=Sys.time()
  mr1=MendelianRandomization::mr_mvivw(object=mvinput)
  hatb1=mr1@Estimate
  stdb1=mr1@StdError
  t1=Sys.time()-tt
  tt=Sys.time()
  mr2=MendelianRandomization::mr_mvegger(object=mvinput,orientate=1)
  hatb2=mr2@Estimate
  stdb2=mr2@StdError.Est
  t2=Sys.time()-tt
  
  tt=Sys.time()
  mr3=MendelianRandomization::mr_mvmedian(mvinput,iterations=100)
  hatb3=mr3@Estimate
  stdb3=mr3@StdError
  t3=Sys.time()-tt
  
  tt=Sys.time()
  tryCatch({mr4=MendelianRandomization::mr_mvlasso(mvinput)},error = function(e) {
    message("An error occurred: ", e)
  })
  hatb4=mr4@Estimate
  stdb4=mr4@StdError
  t4=Sys.time()-tt
  
  tt=Sys.time()
  tryCatch({mr5=MendelianRandomization::mr_mvcML(object=mvinput,rho_mat = Rxy,n=50000,DP=F,num_pert=1)},
           error = function(e) {
             message("An error occurred: ", e)
           })
  hatb5=mr5@Estimate
  stdb5=mr5@StdError
  t5=Sys.time()-tt
  
  tt=Sys.time()
  #tryCatch({mr6=MendelianRandomization::mr_mvcML(object=mvinput,rho_mat=Rxy,n=50000,DP=T,num_pert=100)},
  #error = function(e) {
  #message("An error occurred: ", e)
  #})
  mr6=mr5
  hatb6=mr5@Estimate
  stdb6=mr5@StdError
  t6=Sys.time()-tt
  
  hatb5=mr6@Estimate
  stdb5=mr6@StdError
  t5=Sys.time()-tt
  
  RES11[iter,]=c(hatb1[1],hatb2[1],hatb3[1],hatb4[1],hatb6[1],hatb5[1],hatb0[1])
  RES12[iter,]=c(hatb1[2],hatb2[2],hatb3[2],hatb4[2],hatb6[2],hatb5[2],hatb0[2])
  RES13[iter,]=c(hatb1[3],hatb2[3],hatb3[3],hatb4[3],hatb6[3],hatb5[3],hatb0[3])
  RES14[iter,]=c(hatb1[4],hatb2[4],hatb3[4],hatb4[4],hatb6[4],hatb5[4],hatb0[4])
  RES21[iter,]=c(stdb1[1],stdb2[1],stdb3[1],stdb4[1],stdb6[1],stdb5[1],stdb0[1])
  RES22[iter,]=c(stdb1[2],stdb2[2],stdb3[2],stdb4[2],stdb6[2],stdb5[2],stdb0[2])
  RES23[iter,]=c(stdb1[3],stdb2[3],stdb3[3],stdb4[3],stdb6[3],stdb5[3],stdb0[3])
  RES24[iter,]=c(stdb1[4],stdb2[4],stdb3[4],stdb4[4],stdb6[4],stdb5[4],stdb0[4])
  
  
  TIME[iter,]=c(t1,t2,t3,t4,t5,t0,t0)
  
  print(iter)
  if(iter %% 20 == 0){
    print(colMeans(RES14[1:iter,]))
    print(colMeans(RES24[1:iter,]-colSD(RES14[1:iter,])))
    print(colMeans(TIME[1:iter,]))
    boxplot(RES14[1:iter,])
    remove(G,g_matrix,E,eatX,eatu,X,y)
    #save.image(filename)
  }
  
}
