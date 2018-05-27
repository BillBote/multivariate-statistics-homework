#T1
#这题协方差阵未知需要估计，采用T方检验
#首先处理了txt，使其能被R读取，修改后的文件在文件夹中
rawData<-read.table('MS-4.txt',stringsAsFactors=F)
Data<-as.data.frame(t(rawData))
names(Data)<-c("weight","height")
mat<-cbind(Data[,1],Data[,2])
#上方代码预处理数据得到样本数据阵
mean<-colMeans(mat)#计算样本均值
n1<-length(mat[,1])
A<-t(mat)%*%mat-n1*mean%*%t(mean)#估计协方差阵
mu<-c(63.64,1615.38)
T<-n1*(n1-1)*t(mean-mu)%*%solve(A)%*%(mean-mu)#计算检验统计量
T/(n1-1)/2*(n1-2)#将T方统计量转为F分布统计量
#比较检验统计量和对应分布的分位点
if(T/(n1-1)/2*(n1-2)-qf(0.95,2,n1-2)>0)print("在0.05的显著性水平下拒绝原假设")else print("在0.05的显著性水平下接受原假设")
if(T/(n1-1)/2*(n1-2)-qf(0.99,2,n1-2)>0) print("在0.01的显著性水平下拒绝原假设")else print("在0.01的显著性水平下接受原假设")

#T2
smokeData<-c(168,85,30,17)
smokeP<-smokeData/sum(smokeData)
#估计协方差阵
sigma<-diag(smokeP)-smokeP%*%t(smokeP)
n2<-c(300)
CI<-function(l){
#利用卡方分布的分位数计算置信区间的上下限
low<-t(l)%*%smokeP-sqrt(t(l)%*%sigma%*%l*qchisq(0.95,3)/n2)
high<-t(l)%*%smokeP+sqrt(t(l)%*%sigma%*%l*qchisq(0.95,3)/n2)
cat("95%的联立置信区间是[",low,",",high,"].")
}
CI(c(1,2,1,2))#例子
