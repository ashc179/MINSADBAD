####--------------------------------------------------------------------####
####-------------------------Regresion MINSADBAD------------------------####
####--------------------------------------------------------------------####

n=40        #Número de observaciones

B0=3
B1=4       #Parametros
B2=5

X1=runif(n,1,4)
X2=runif(n,2,5)
ei=rnorm(1000,0,3)
e=matrix(ei,n,1000)

ef=c(matrix(ei,n,1000),(0.05*e/3)+(1-0.05)*(e),(0.10*e/3)+(1-0.10)*(e),(0.15*e/3)+(1-0.15)*(e)
,(0.20*e/3)+(1-0.20)*(e),(0.2)^(-1)*(exp(0.2*e/3)-1),0.4^(-1)*(exp(0.4*e/3)-1)
,0.6^(-1)*(exp(0.6*e/3)-1),0.8^(-1)*(exp(0.8*e/3)-1),1^(-1)*(exp(1*e/3)-1)
,0.2*exp((0.2*(e/3)^2)/2),0.4*exp((0.4*(e/3)^2)/2),0.6*exp((0.6*(e/3)^2)/2)
,0.8*exp((0.8*(e/3)^2)/2),1*exp((1*(e/3)^2)/2))

err=matrix(ef,n,length(ef)/n)

Yi=rep(0,length(ef))
Yi=matrix(Yi,n,length(ef)/n)

for (i in 1:length(err[,1])){
for (j in 1:length(err[1,])){
Yi[i,j]=B0+B1*X1[i]+B2*X2[i]+err[i,j]}}
Yi

Y=t(t(apply(Yi,1,median)))
Y

##------------------------##
##--------Graficas--------##
##------------------------##

hist(Y,freq=F,xlim=c(15,50),ylim=c(0,0.10))
x_=mean(Y)
ds=sd(Y)
curve(dnorm(x,x_,ds),add=T,col="dodgerblue2")

#install.packages("lmtest")
library(lmtest)
datos=data.frame(X1, X2, Y)
plot(datos)
X=model.matrix( ~ datos$X1 + datos$X2)
fit1=lm(formula = datos$Y ~ datos$X1 + datos$X2, data = datos)
summary(fit1)

source("C:/Tareas/Funcionesauxiliares1.txt", encoding = "Latin1")
#source("C:\\Users\\cnaber\\Trabalho\\windows\\Unicamp\\Disciplinas\\2_semestre_2016\\ME 613\\Programas\\diag_norm.r")
envelnorm(fit1)
#source("http://www.poleto.com/funcoes/diag.norm.txt")
diagnorm(fit1)

##------------------##
##-----Benites------##
##------------------##

#install.packages("xtable")
library(xtable)
#install.packages("lqr")
library(lqr)

res.aju.mod.lm=function(fit.model,gama,digitos=4)
{
  mX=model.matrix(fit1)
  sumres=summary(fit1)
  n=nrow(mX)
  p=ncol(mX)
  quantilt=qt(0.5*(1+gama),df=n-p)
  v.beta=fit.model$coefficients
  ep.beta=sqrt(diag(vcov(fit.model)))
  result=cbind(sumres$coefficients[,1:2],v.beta-quantilt*ep.beta,v.beta+quantilt*ep.beta,sumres$coefficients[,3],sumres$coefficients[,4])
  colnames(result)=c("Est.","EP","LIIC","LSIC","Estat. t","p-valor")
  print(xtable(result,digits=digitos))
  resultl=list(result=result)
  return(resultl)
}
res.aju.mod.lm(fit1,0.05,4)
ppp=par(mfrow=c(2,2))
plot(fit1)
par(ppp)
r=resid(fit1)
shapiro.test(r)
p=0.5
fit2=best.lqr(datos$Y ~ datos$X1 + datos$X2,p=0.5,precision=10^-6, criterion="AIC")
##------------------------##
##--------Matriz H--------##
##------------------------##

H=cbind(matrix(1,n-1),-diag(n-1))
for (i in 1:(n-2)){
  h = cbind(matrix(0,n-(i+1),i),matrix(1,n-(i+1)),-diag(n-(i+1)))
  H=rbind(H,h)
}
H

HY=H%*%Y
HY

##------------------------##
##--------Sujeto a--------##
##------------------------##

p=2               #Número de variables independientes
r=(n*(n-1))/2     #Restricciones
r

X=cbind(X1,X2)
X

S1=cbind(X,diag(n),matrix(0,n,r),matrix(0,n,r))
S2=cbind(-X,diag(n),matrix(0,n,r),matrix(0,n,r))
S3=cbind(matrix(0,r,p),H,diag(r),-diag(r))
S=rbind(S1,S2,S3)
S

b11=t(t(Y[,1]))
b21=t(t(-Y[,1]))
b31=matrix(0,r)
b=rbind(b11,b21,b31)
b

fun.obj=rbind(0,0,matrix(0,n),matrix(1,r),matrix(1,r))
fun.obj

lb=rep(0,length(fun.obj))
lb

##----------------------------------------------------------##
##---Método Simplex para problemas de Programación Lineal---##
##----------------------------------------------------------##

#install.packages("boot")
library(boot)

a=c(fun.obj)
A1=c(rep(0,length(fun.obj)))
b1=c(0)
A2=rbind(S1,S2)
b2=c(b11,b21)
A3=S3
b3=c(b31)

x=simplex(a=a,A1,b1,A2,b2,A3,b3,maxi=FALSE)
x

##----------------------------------##

B1est=x$soln[1]
B1est
B2est=x$soln[2]
B2est

##------------------------##
##----Estimación de ß0----##
##------------------------##

datos=data.frame(X1, X2, Y)
X=model.matrix( ~ datos$X1 + datos$X2)
X
B0esti=solve(t(X)%*%X)%*%t(X)%*%Y
B0est=B0esti[1,1]
B0est

Best=rbind(B0est,B1est,B2est)
Best

##-------------------------##
##----B1 y B2 estimados----##
##-------------------------##

a=seq(1,15000,by=1000)
b=seq(1000,15000,by=1000)

Yy=matrix(0,n,15)
for(k in 1:15){ 
Yy[,k]=t(t(apply(Yi[,a[k]:b[k]],1,median)))
}
Yy

a=c(fun.obj)
A1=c(rep(0,length(fun.obj)))
b1=c(0)
A2=rbind(S1,S2)
A3=S3
b31=matrix(0,r)
b3=c(b31)

b11=Yy
b21=-Yy
b31=matrix(0,r,15)
b=rbind(b11,b21,b31)

bb=rbind(b11,b21)

B1est=rep(0,15)
B2est=rep(0,15)
for(i in 1:15){
x=simplex(a=a,A1,b1,A2,bb[,i],A3,b3,maxi=FALSE)
B1est[i]=x$soln[1]
B2est[i]=x$soln[2]
}
t(t(B1est))
t(t(B2est))

##------------------------##
##----Estimación de ß0----##
##------------------------##

X=model.matrix( ~ datos$X1 + datos$X2)
X

B0est=rep(0,15)
for(i in 1:15){
B0es=solve(t(X)%*%X)%*%t(X)%*%t(t(Yy[,i]))
B0est[i]=B0es[1,1]
}
t(t(B0est))

##----------------------------------##
##-----Criterios de Comparación-----##
##----------------------------------##

##----------##
##----ECM---##
##----------##

B0=rep(B0,15)
B1=rep(B1,15)
B2=rep(B2,15)
B=rbind(B0,B1,B2)
B

Best=rbind(B0est,B1est,B2est)
Best

ECM=t(t(apply((Best-B)^2,2,mean)))
ECM

ECM0=t(t(apply((t(B0est)-B0)^2,2,mean)))
ECM0

ECM1=t(t(apply((t(B1est)-B1)^2,2,mean)))
ECM1

ECM2=t(t(apply((t(B2est)-B2)^2,2,mean)))
ECM2

##----------##
##----Ea----##
##----------##

Ea=t(t(apply(abs(B-Best),2,mean)))
Ea

Ea0=t(t(apply(abs(B0-t(B0est)),2,mean)))
Ea0
 
Ea1=t(t(apply(abs(B1-t(B1est)),2,mean)))
Ea1

Ea2=t(t(apply(abs(B2-t(B2est)),2,mean)))
Ea2

##----------##
##----Er----##
##----------##

Er=t(t(apply(abs(B-Best)/abs(B),2,mean)))
Er

Er0=t(t(apply(abs(B0-t(B0est))/abs(B0),2,mean)))
Er0

Er1=t(t(apply(abs(B1-t(B1est))/abs(B1),2,mean)))
Er1

Er2=t(t(apply(abs(B2-t(B2est))/abs(B2),2,mean)))
Er2

##----------------------------------------------##
##------Grafica - Criterios de comparación------##
##----------------------------------------------##

#install.packages("reshape",dependencies=TRUE)
#install.packages("ggplot2")

##--------------------------##
##-------Contaminadas-------##
##--------------------------##

x=c(0.05,0.10,0.15,0.20)
y1=c(ECM[2:5,])
y2=c(Ea[2:5,])
y3=c(Er[2:5,])

df=data.frame(X=x,ECM=y1,Ea=y2,Er=y3)
df

library(reshape)
db=melt(df,id="X")
names(db)=c("X","met","Y")
db

tab=rep("Contaminadas",12)
db1=cbind(db,tab)
db1

##--------------------------##
##-------g de Tukey---------##
##--------------------------##

x=c(0.2,0.4,0.6,0.8,1)
y1=c(ECM[6:10,])
y2=c(Ea[6:10,])
y3=c(Er[6:10,])

df=data.frame(X=x,ECM=y1,Ea=y2,Er=y3)
df

library(reshape)
db=melt(df,id="X")
names(db)=c("X","met","Y")
db

tab=rep("g-Tukey",15)
db2=cbind(db,tab)
db2

##--------------------------##
##--------h de Tukey--------##
##--------------------------##
x=c(0.2,0.4,0.6,0.8,1)
y1=c(ECM[11:15,])
y2=c(Ea[11:15,])
y3=c(Er[11:15,])

df=data.frame(X=x,ECM=y1,Ea=y2,Er=y3)
df

library(reshape)
db=melt(df,id="X")
names(db)=c("X","met","Y")
db

tab=rep("h-Tukey",15)
db3=cbind(db,tab)
db3

##-----Pegando bases--------##

db=rbind(db1,db2,db3)
names(db)

##--------Grafico-----------##

library(ggplot2)
p=ggplot(db,aes(x=X,y=Y,colour=met,group=met)) + facet_wrap(~tab,scales = "free_x")
pg=p + geom_point() + geom_line() + labs(colour="Criterios de Comparación",x="Distribuciones", y="Criterios de Comparación", title="Regresión Minsadbad") +scale_y_log10()
pg

##----------##
##----DC----##
##----------##

Dc=cooks.distance(fit1)
D=as.vector(Dc)
D
p=2

cota=4/(n-p-2)   #Para identificar los puntos influyentes en el modelo 

significativas=Dc>cota
significativas   
##Ninguna observación supera el valor referencia (cota), 
##ninguna observación es potencialmente influyente en el modelo

plot(fit1,which=4,cook.levels=cota,las=1)
abline(h=cota,lty="dashed",col="dodgerblue2")

par(mfrow=c(2,2))   #Residuales del fit1
plot(fit1, col='deepskyblue4', pch=19)

##----------##
##----DA----##          
##----------##

X=model.matrix( ~ X1 + X2)
X
H=X%*%solve(t(X)%*%X)%*%t(X) 
h=diag(H)
h
lms=summary(fit1)
s=lms$sigma
r=resid(lms)
ti=r/(s*(1-h)^(1/2))
tsi=ti*((n-p-1)/(n-p-ti^2))^(1/2)

Ai=((((n-p)/p)*(h/(1-h)))^(1/2))*abs(tsi)
A=as.vector(Ai)
A

cota=2*(p/(n-p))^(1/2)   #Para identificar los puntos influyentes en el modelo 

significativas=Ai>cota
significativas   
##()()() son observación que superan el valor referencia (cota), 
##()()() son observaciones potencialmente influyente en el modelo

plot(A,type="h",main= "Atkinson distance",ylab="Atkinson distance",xlab="Obs. number")
abline(h=cota, lty="dashed", col="dodgerblue2")

#########################################################
###Histograma residuos estudentizados y estandarizados###
#########################################################

hist(ti,freq=F,xlim=c(-3,3),ylim=c(0,0.8))
x_=mean(ti)
ds=sd(ti)
curve(dnorm(x,x_,ds),add=T,col="dodgerblue2")
qqnorm(ti)
qqline(ti)

hist(tsi,freq=F,xlim=c(-4,4),ylim=c(0,0.8))
x_=mean(tsi)
ds=sd(tsi)
curve(dnorm(x,x_,ds),add=T,col="dodgerblue2")
qqnorm(tsi)
qqline(tsi)

###################

datos=data.frame(B0est,B1est,B2est,ECM,ECM0,ECM1,ECM2,Ea,Ea0,Ea1,Ea2,Er,Er0,Er1,Er2)
distancias=data.frame(D,A)
#install.packages("openxlsx")
library(openxlsx)
write.xlsx(datos,"Resultados40.xlsx")
write.xlsx(distancias,"Resultados40.xlsx")
read.xlsx(datos)
read.xlsx(distancias)

