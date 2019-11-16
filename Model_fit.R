#*****GAMMAR_fit.R �ù���ӻ��ģ���������*******
library(foreign)
library(splines)
library(mgcv)
library("tseries")

testdata<-read.csv("F:\\Jinger\\test729.csv",header=T)
testdata=testdata[1:54,]
path="F:\\Jinger\\R_programs"
setwd(path) 
source('GAMAR.R')
##*******���ģ��********
glmbb <- gam((case)~ns(temperature,2)+ns(humidity,3)+ns(precipitation,3)+ns(index,3)+year+month,data=testdata,method="GCV.Cp",family=poisson(link = "log"))
result <- gamAR(case~ns(temperature,2)+ns(humidity,3)+ns(precipitation,3)+ns(index,3)+year+month,data=testdata,p.ar=1,starts=glmbb$coef);
write.csv(result$coefficients,file="F:/Jinger/�����ļ�/�Ǹ���/test/Estimate.csv",row.names=F,quote=F);
#write.csv(result$X,file="F:/Jinger/�����ļ�/�Ǹ���/test/s1.csv",row.names=F,quote=F);
#write.csv(result$fitted.values,file="F:/Jinger/�����ļ�/�Ǹ���/test/fitted_data.csv",row.names=F,quote=F);

casedata=as.matrix(testdata$case)
total_num=nrow(casedata)

#*******gamAR�������********
fit_data=as.matrix(result$fit)
row_num=nrow(fit_data)

#********�������������ʵ�����ݵĶԱ�ͼ**********
data1=cbind(casedata[total_num-row_num+1:total_num],fit_data[,1])
colnames(data1) <- c("truedata","fitteddata")
fit.int=lm(data1[,1] ~ data1[,2])
dev.new();
plot(data1[,1],data1[,2],xlab="Observation",ylab="fitted data")
abline(fit.int)

#*******�����ֵ��ʵ��ֵ��ʱ���ͼ********
dev.new()
x=seq(as.Date("2011/4/1"), as.Date("2015/6/1"), by="1 months")
plot(x,log10(data1[1:51,1]),type="l",col="blue",lty=4,xlab="Time",ylab="Case")
#lines(x,log10(data1[1:51,1]),type="l",lty=4,col="blue",xlab="Time",ylab="Case")
lines(x,log10(data1[1:51,2]),lty=2,col="red",xlab="Time",ylab="Case")
text.legend=c("Observation","Fitted data")
col2<-c("blue","red")
legend("topleft",legend=text.legend,col=c(col,col2),lty=c(4,2),bty="n",horiz=FALSE,trace=TRUE,plot=TRUE)

#*******Peason�в������غ�ƫ���ͼ*********
dev.new()
attach(mtcars)
par(mfrow=c(2,1))
Pearson_residuals=result$Pearson.residuals
acf(Pearson_residuals,lag.max=4) #�����
pacf(Pearson_residuals,lag.max=4) #ƫ���ͼ

#*********��Pearson�в���Ԥ��ֵɢ��ͼ********
dev.new();
plot(result$fit,result$Pearson.residuals,xlab="predicted data",ylab="Pearson residuals")

#*********��Pearson�в���ʱ��ɢ��ͼ**********
dev.new();
x=seq(as.Date("2011/3/1"), as.Date("2015/6/1"), by="1 months")
plot(x,result$Pearson.residuals[1:52],xlab="Time",ylab="Pearson residuals",ylim=c(-15,15))