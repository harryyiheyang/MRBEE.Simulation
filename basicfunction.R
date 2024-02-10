colSD=function(A){
B=A
p=length(A[1,])
for(i in 1:p){
B[,i]=sd(A[,i])  
}
return(B)
}

vecnorm=function(a){
sqrt(sum(a^2))
}

CScov=function(p,rho){
A=diag(p)*(1-rho)+matrix(rho,p,p)
return(A)
}

standard=function(A){
b=apply(A,2,sd)
a=colMeans(A)
A=t(t(A)-a)
A=t(t(A)/b)
return(A)
}

Uerrorvar=function(Sbb,hxvec){
p=length(hxvec)
Vb=diag(Sbb)
D=Vb/hxvec
return(D)
}

verrorvar=function(Sbb,Ruv,Du,theta,hy){
p=length(theta)
Suu=t(Ruv[1:p,1:p]*sqrt(Du))*sqrt(Du);
ruv=Ruv[1:p,p+1]
svv=mvvarv(Vbb=Sbb,Vuu=Suu,Ruv=ruv,theta=theta,hy=hy)  
return(svv)
}

CovXy=function(Sbb,Suv,Rover,theta){
p=length(theta)
Suu=Suv[1:p,1:p]
Svv=Suv[p+1,p+1]
Suv=Suv[1:p,p+1]
Sxy=diag(p+1)
Sxy[1:p,1:p]=Suu+Sbb
Sxy[1:p,p+1]=Sxy[p+1,1:p]=(Suu+Sbb)%*%theta+Suv
Sxy[1+p,1+p]=t(theta)%*%(Suu+Sbb)%*%theta+2*sum(theta*Suv)+Svv
Vxy=Sxy*Rover 
return(Vxy)
}

rangeinput=function(theta,frac=0.7){
A=cbind(theta*frac,theta/frac)
B=apply(A,1,sort)  
return(B)
}

biggwas=function(x,G){
x=as.vector(x)
ux=mean(x)
vx=var(x);vx=as.numeric(vx)
ug=colMeans(G)
G=t(t(G)-ug)
vg=colSums(G^2)
b=(t(G)%*%(x-ux))/vg
sdb=(vx-b^2*vg/length(x))/length(x)
A=list()
A$est=as.vector(b)
A$std=as.vector(sqrt(sdb))
return(A)
}

mbiggwas=function(X,G1,N){
p=dim(X)[2]
m=dim(G1)[2]
B=matrix(0,m,p)
SDB=B
for(i in 1:p){
x=X[N[,i],i]
G=G1[N[,i],]
x=as.vector(x)
ux=mean(x)
vx=var(x);vx=as.numeric(vx)
ug=colMeans(G)
G=t(t(G)-ug)
vg=colSums(G^2)
b=matrixVectorMultiply(t(G),x-ux)/vg
b=as.vector(b)
sdb=(vx-b^2*vg/length(x))/length(x)
B[,i]=b
SDB[,i]=sqrt(as.vector(sdb))
}
A=list()
A$est=B
A$std=SDB
return(A)
}


validadj <- function(vector1, vector2, tau) {
diff <- length(vector2) / length(vector1)

if (diff < tau) {
missing_indices <- setdiff(1:length(vector1), vector2)
sorted_missing_indices <- missing_indices[order(vector1[missing_indices])]
num_to_add <- ceiling(tau * length(vector1)) - length(vector2)
vector2 <- c(vector2, sorted_missing_indices[1:num_to_add])
}

return(vector2)
}
