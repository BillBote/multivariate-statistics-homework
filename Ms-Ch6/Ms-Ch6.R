#T1
rawdata1<-read.table("MS-6.txt")
rawdata1<-rawdata1[,-1]
names(rawdata1)<-c("x1","x2","x3","x4","type")
mat_1<-rawdata1[which(rawdata1$type==1),]
mat_2<-rawdata1[which(rawdata1$type==2),]
mat_1<-mat_1[,-5]
mat_2<-mat_2[,-5]
mu_1<-colMeans(mat_1)
mu_2<-colMeans(mat_2)
mat_1<-cbind(mat_1[,1],mat_1[,2],mat_1[,3],mat_1[,4])
mat_2<-cbind(mat_2[,1],mat_2[,2],mat_2[,3],mat_2[,4])
n1<-length(mat_1[,1])
n2<-length(mat_2[,1])
sigma_1<-(t(mat_1)%*%mat_1-n1*mu_1%*%t(mu_1))/(n1-1)
sigma_2<-(t(mat_2)%*%mat_2-n2*mu_2%*%t(mu_2))/(n2-1)
#距离判别
DDiscrimnation<-function(x,mu1,mu2,sigma1,sigma2){
  if(length(mu1)!=length(mu2))
    print("Two mean must have same length!")
  if(length(x)!=length(mu1))
    print("x must have same length with means!")
  distance<-t(x-mu1)%*%solve(sigma1)%*%(x-mu1)-t(x-mu2)%*%solve(sigma2)%*%(x-mu2)
  if(distance<=0)
    class<-1
  else
    class<-2
  return(class)
}
#本题的一个例子
x<-c(169,287,162,214)
class<-DDiscrimnation(x,mu_1,mu_2,sigma_1,sigma_2)
cat("x belongs to class",class,"!")

#bayes判别
bayesDiscrimnation<-function(x,mu1,mu2,sigma1,sigma2){
  if(length(mu1)!=length(mu2))
    print("Two mean must have same length!")
  if(length(mu1)!=length(x))
    print("x must have same length with means!")
  W<--1/2*t(x)%*%(solve(sigma1)-solve(sigma2))%*%x+(t(mu1)%*%solve(sigma1)-t(mu2)%*%solve(sigma2))%*%x
  k<-log((1/3)/(2/3))+1/2*log(det(sigma1)/det(sigma2))+1/2*(t(mu1)%*%solve(sigma1)%*%mu1-t(mu2)%*%solve(sigma2)%*%mu2)
  if(W>=k)
    class<-1
  else
    class<-2
  return(class)
}
#一个例子
x<-c(169,287,162,214)
class<-bayesDiscrimnation(x,mu_1,mu_2,sigma_1,sigma_2)
cat("x belongs to class",class,"!")

#比较明显误判率
sum1<-0
sum2<-0
#距离判别误差率
class1<-c()
for(i in 1:length(rawdata1[,1])){
  class1[i]<-DDiscrimnation(as.numeric(rawdata1[i,][-5]),mu_1,mu_2,sigma_1,sigma_2)
  if(class1[i]!=as.numeric(rawdata1[i,][5]))
    sum1<-sum1+1
  p1<-sum1/length(rawdata1[,1])
}
as.data.frame(class1)#查看具体分类情况
#bayes判别误差率
class2<-c()
for(i in 1:length(rawdata1[,1])){
  class2[i]<-bayesDiscrimnation(as.numeric(rawdata1[i,][-5]),mu_1,mu_2,sigma_1,sigma_2)
  if(class2[i]!=as.numeric(rawdata1[i,][5]))
    sum2<-sum2+1
  p2<-sum2/length(rawdata1[,1])
}
as.data.frame(class2)#查看具体分类情况
#比较
p1
p2#查看明显误分率

#T2
rawdata2<-read.table("MS-6(2).txt")
rawdata2<-as.data.frame(t(as.matrix(rawdata2)))
names(rawdata2)<-c("x1","x2","type")
c<-cbind(c(0,10,50),c(500,0,200),c(100,50,0))
p<-c(0.05,0.6,0.35)
mat1<-as.matrix(rawdata2[which(rawdata2$type==1),][,-3])
mat2<-as.matrix(rawdata2[which(rawdata2$type==2),][,-3])
mat3<-as.matrix(rawdata2[which(rawdata2$type==3),][,-3])
mu1<-colMeans(mat1)
mu2<-colMeans(mat2)
mu3<-colMeans(mat3)
mu<-cbind(mu1,mu2,mu3)
n1<-length(mat1[,1])
n2<-length(mat2[,1])
n3<-length(mat3[,1])
sigma<-(t(mat1)%*%mat1-n1*mu1%*%t(mu1)+t(mat2)%*%mat2-n2*mu2%*%t(mu2)+t(mat3)%*%mat3-n3*mu3%*%t(mu3))/(n1+n2+n3-3)

#bayes判别函数
bayes<-function(x,p,c,mu,sigma){
  f<-c()
  h<-c()
  for(i in 1:3)
    f[i]<-1/(2*pi)*1/sqrt(det(sigma))*exp(-1/2*t(x-mu[,i])%*%solve(sigma)%*%(x-mu[,i]))
  for(i in 1:3){
    h[i]<-0
    for(j in 1:3)
      h[i]<-h[i]+p[j]*c[i,j]*f[j]
  }
  class<-which.min(h)
  return(class)
}

#fisher判别
B<-0
mean<-rowMeans(mu)
for(i in 1:3){
  B<-B+(mu[,i]-mean)%*%t(mu[,i]-mean)
}
library(pracma)
sqrtB<-sqrtm(B)$B
lambda<-eigen(B)
vector<-lambda$vectors
alpha1<-vector[,1]
alpha2<-vector[,2]
l1<-solve(sqrtB)%*%alpha1
l2<-solve(sqrtB)%*%alpha2
y<-matrix(nrow=2,ncol=3)
for(i in 1:3){
  y[,i]<-c(t(l1)%*%mu[,i],t(l2)%*%mu[,i])
}
#计算欧式距离
distance2D<-function(x,y){ 
  sum<-0
  for(i in 1:length(x)){
    sum<-sum+(x[i]-y[i])^2
  }
  return(sum)
}
#fisher判别函数
fisher<-function(x){
  yx<-c(t(l1)%*%x,t(l2)%*%x)
  Dis<-c()
  for(i in 1:3){
    Dis[i]<-distance2D(yx,y[,i])
  }
  class<-which.min(Dis)
  return(class)
}

#bayes的明显误分率
Class1<-c()
SUM1<-0
for(i in 1:length(rawdata2[,1])){
  Class1[i]<-bayes(as.numeric(rawdata2[i,][-3]),p,c,mu,sigma)
  if(Class1[i]!=as.numeric(rawdata2[i,][3]))
    SUM1<-SUM1+1
  P1<-SUM1/length(rawdata2[,1])
}
P1
Class1

#fisher误分率
Class2<-c()
SUM2<-0
for(i in 1:length(rawdata2[,1])){
  Class2[i]<-fisher(as.numeric(rawdata2[i,][-3]))
  if(Class2[i]!=as.numeric(rawdata2[i,][3]))
    SUM2<-SUM2+1
  P2<-SUM2/length(rawdata2[,1])
}
P2
#检测样本
x<-rbind(c(195,123),c(211,122),c(187,123),c(192,109))
Class3B<-c()
Class3F<-c()
for(i in 1:length(x[,1])){
  Class3B[i]<-bayes(x[i,],p,c,mu,sigma)
  Class3F[i]<-fisher(x[i,])
}
Class3B
Class3F#查看两种判别的分类情况
