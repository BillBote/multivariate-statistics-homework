data<-read.table("MS-Ch3.txt",head=T)#我重写了txt文件用来保存数据，重写的文件放在文件夹里
n<-27
mat_1<-cbind(data$X1971.1,data$X1971.2,data$X1971.3,data$X1971.4,data$X1971.5,data$X1971.6)
mat_2<-cbind(data$X1972.1,data$X1972.2,data$X1972.3,data$X1972.4,data$X1972.5,data$X1972.6)#这里可以用as.matrix来转为矩阵
mu_1<-colMeans(mat_1)
mu_2<-colMeans(mat_2)
mu_1
mu_2

A_1<-t(mat_1)%*%mat_1-n*mu_1%*%t(mu_1)
A_2<-t(mat_2)%*%mat_2-n*mu_2%*%t(mu_2)
sigma_1<-A_1/n
sigma_2<-A_2/n
sigma_1
sigma_2

cor_1<-cor(mat_1)
cor_2<-cor(mat_2)
cor_1
cor_2
