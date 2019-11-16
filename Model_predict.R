library(foreign)
library(splines)
library(mgcv)
library("tseries")

testdata<-read.csv("F:\\Jinger\\test729.csv",header=T)
Estimatecoeff<-read.csv("F:\\Jinger\\静儿文件\\登革热\\test\\Estimate.csv",header=T)

path="F:\\Jinger\\R_programs"
setwd(path) 
source('GAMAR.R')

#*******求解模型********
glmbb <- gam((case)~ns(temperature,2)+ns(humidity,3)+ns(precipitation,3)+ns(index,3)+year+month,data=testdata,method="GCV.Cp",family=poisson(link = "log"))
result <- gamAR(case~ns(temperature,2)+ns(humidity,3)+ns(precipitation,3)+ns(index,3)+year+month,data=testdata,p.ar=1,starts=glmbb$coef);
write.csv(result$X,file="F:/Jinger/静儿文件/登革热/test/s1.csv",row.names=F,quote=F);


coeff=as.matrix(Estimatecoeff[1:12,1]);
sdata=as.matrix(sdata);
nsdata=sdata%*%coeff
ARcoeff=as.matrix(Estimatecoeff[13:15,1]);
data=as.matrix(truedata$case);
yy=matrix(1:18,nrow=6,ncol=3);
bb=matrix(1:18,nrow=6,ncol=3);

for(i in 1:6)
{     
	if(data[54+i-1,1]>0.5)
          yy[i,1]=data[54+i-1,1] else
           yy[i,1]=0.5;
           bb[i,1]= nsdata[54+i-1,1];  
      if(data[54+i-2,1]>0.5)
          yy[i,2]=data[54+i-2,1] else
           yy[i,2]=0.5;
           bb[i,2]= nsdata[54+i-2,1]; 
      if(data[54+i-3,1]>0.5)
          yy[i,3]=data[54+i-3,1] else
           yy[i,3]=0.5;
          bb[i,3]= nsdata[54+i-3,1]; 
     
}
cc=log(yy)-bb;
AR=cc%*%ARcoeff;
#w=rep(1:6);
predict=nsdata[55:60]+AR;
predict=exp(predict)


write.csv(predict,file="F:/Jinger/静儿文件/登革热/test/predicteddata.csv",row.names=F,quote=F)



