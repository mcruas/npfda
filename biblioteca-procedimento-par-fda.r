library(fda)

logprecav = CanadianWeather$dailyAv[dayOfYearShifted, ,'log10precip']
plot.ts(logprecav[,1], col = 1)


vetor=c(1,2,-3,3,5,10,3,5,4,3,2,5)
dayrange = c(0,365)
basis.spline.365 = create.bspline.basis(c(0,12),length(vetor))

coef=eval.basis(c(0,0.2,0.5),basis.spline.365)



sinefd=fd(vetor,basis.spline.365)
plot(sinefd)

x <- function x^2


plot(basis.spline.365)

sinefd$basis

plot(tst)



daybasis = create.fourier.basis(dayrange,365)



Lcoef = c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd = vec2Lfd(Lcoef,dayrange)


basis13 = create.bspline.basis(c(0,1),13)
tvec = seq(0,1,len=5)
sinecoef=sin(2*pi*tvec)
sinefd=fd(sinecoef,basis13,fdnames=list("a","b","c"))
op = par(cex = 1.2)
plot(sinefd,lwd=2)
points(tvec,sinecoef,lwd=2)
par(op)

testebase = create.bspline.basis(rangeval=c(0,6),6)  
blau = eval.basis(c(1,2,5,1,2,10),basisobj=testebase)
tst = fd(blau,testebase)
plot(tst)
objfd <- fd(coef=blau,basisobj=testebase)

# Cria uma base para bsplines
# o segundo argumento é o número de funcoes base = ordem (4, normalmente) + número de nós interiores
basisobj = create.bspline.basis(rangeval=c(0,1),nbasis=(2 + dim(y.past.test)[2]))  

#fd cria um objeto funcional
tstfn <- fd(y.past.test[,1],basisobj)

plot(basisobj)

plot(tstfn)

plot(y,type="l")
dados = parte.dados.ettj[,1]

# primeiro argumento são os dados, segundo a base de funções
spline = eval.basis(evalarg=,basisobj=basisobj)


estimacao.linear <- function(base_nefasta,maturidade,s,percentual.testar,truncar,obs.por.curva,todas.maturidades)
{
  if(passo < 1){
    stop(paste("Você forneceu um passo inválido"))
  }
  
  z <- ts(data=base_nefasta[truncar[1]:truncar[2],maturidade+1])
  tamanho.serie.temporal  = length(z)
  
  learning <- c(1,trunc(percentual.testar*tamanho.serie.temporal)-obs.por.curva-s+1)
  testing <- c(learning[2],tamanho.serie.temporal-s-obs.por.curva+1)
  
  y.past.learn <- NULL
  y.past.test <- NULL
  y.futur.s <- NULL
  y.test.true.values <- NULL
  
  i = learning[2] - s - obs.por.curva + 1 # o índice i representa a primeira observação das curvas que serão armazenadas
  #laço vai preenchendo as curvas de tras para frente
  while (i >= learning[1]) 
  {
    y.past.learn <- cbind(z[i:(i+obs.por.curva-1)],y.past.learn)  
    y.futur.s <- c(z[i+obs.por.curva-1+s],y.futur.s)
    i = i - passo 
  }
  y.past.learn <- t(y.past.learn)
  
  # o índice j é índice da matriz y, enquanto i é índice da série temporal z
  for (i in testing[1]:testing[2])
  {
    y.past.test <- cbind(y.past.test,z[i:(i+obs.por.curva-1)]) 
    y.test.true.values <- c(y.test.true.values,z[i+obs.por.curva-1+s])
  }
  y.past.test <- t(y.past.test)
  
  ##
  y.test.true.values.2 <- y.test.true.values
  ##
  
  media.mes <- rowMeans(y.past.learn, na.rm = TRUE)
  media.test <- rowMeans(y.past.test, na.rm = TRUE)
  
  y.past.learn.media <- y.past.learn - media.mes
  y.futur.s.media <- y.futur.s - media.mes
  y.past.test.media <- y.past.test - media.test
  
  result.pred.s.media <- funopare.knn.lcv(y.futur.s.media,y.past.learn.media,y.past.test.media,n.componentes.principais,kind.of.kernel="quadratic",semimetric="pca")
  result.pred.s.pca <- result.pred.s.media$Predicted.values + media.test
  
  result.pred.s.media <- funopare.knn.lcv(y.futur.s.media,y.past.learn.media,y.past.test.media,grau.derivada,nknot=6, c(0,20),kind.of.kernel="quadratic",semimetric="deriv")
  result.pred.s.deriv <- result.pred.s.media$Predicted.values + media.test
  
  
  mse.ar = (y.past.test[,obs.por.curva]-y.test.true.values)^2
  mse.pca = (result.pred.s.pca - y.test.true.values)^2
  mse.deriv = (result.pred.s.deriv - y.test.true.values)^2
  mse.average.pca.ar = ((1/2)*(result.pred.s.pca + y.past.test[,obs.por.curva])-y.test.true.values)^2
  
  media.mse.ar = mean(mse.ar)
  media.mse.pca = mean(mse.pca)
  media.mse.deriv = mean(mse.deriv)
  media.mse.average.pca.ar = mean(mse.average.pca.ar)
  
  data=date()
  linha=paste(data,media.mse.average.pca.ar,media.mse.ar,media.mse.pca,media.mse.deriv,maturidade,obs.por.curva,s,n.componentes.principais,grau.derivada,truncar[1],truncar[2],passo,sep=',')
  
}