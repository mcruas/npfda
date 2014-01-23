
# Funções --------------------------------------------------------------


LocArq <- function (dir, intervalo, hor, metodo, maturidade) {
  return(paste0(dir,intervalo[1],":",intervalo[2],"-",hor,"-",metodo,"-",maturidade))
}

nome.vetor <- function(object) {
  return(paste0(object[1],":",object[2]))
}


Plot.list <- function(object) {
  tmp <- as.data.frame(object)
  Multiplot.ts(tmp)
}

EQM <- function(d1,d2) {
  mean((d1-d2)^2)
}


# Parâmetros para a estimação ---------------------------------------------

maturidade <- 5 # maturidade do titulo
obs.por.curva = 50 #quantas observações por curva
percentual.testar = 0.2 #quantos % da amostra usar para treinamento. o resto será para teste
#quais.obs = c(1,20) # quais observações para mostrar no gráfico
s <- 128 # horizonte de planejamento
n.componentes.principais <- 5 # quais componentes principais testar
grau.derivada <- 0 # qual grau da derivada testar
truncar <- c(2500,5500) #um número que mostra quais os limites que são analisados da serie temporal: deve variar entre 1 e 6537
passo <- 50
todas.maturidades <- 1:20


## Diferentes opções de parâmetros para testes
# vetor.truncar <- list(c)
# vetor.grau.derivada = c(0,1) #quais derivadas testar 
# vetor.n.componentes.principais = c(7) # quais componentes principais testar
# vetor.s = c(128)
# todas.maturidades <- 1:20
# vetor.obs.por.curva = c(30,40,50,60,100,120)
# vetor.n.vizinhos = c(2:15*2)
# vetor.maturidade = c(1)
# vetor.passo = c(15,20,30,40,60,80,90)

## Testes totais
vetor.truncar <- list(c(1,1000)) # ,c(201,1200),c(401,1400),c(601,1600)
vetor.grau.derivada = c(0) #quais derivadas testar 
vetor.n.componentes.principais = c(7) # quais componentes principais testar
vetor.s = c(10,32)
vetor.obs.por.curva = c(50, 60) # 30,40
vetor.n.vizinhos = c(1:10*2)
vetor.maturidade = c(1:4,6,8,12,17)
vetor.passo = c(5,15,20,30)


intervalos <- list(c(1,1500),c(501,2000),c(1001,2500),c(1,2500))
vetor.grau.derivada = c(0,1,2) #quais derivadas testar 
vetor.n.componentes.principais = c(5,6,7) # quais componentes principais testar
vetor.s = c(21,125,250)
vetor.obs.por.curva = c( 50) # 30,40, 50, 60,
vetor.n.vizinhos = c(1:10*2)
vetor.maturidade = c(1:4,6,8,12,17)
vetor.passo = c(5,15,20,30)


# para teste
intervalos <- list(c(1,100),c(51,151))
horizonte <- 8

# Código principal --------------------------------------------------------

#library(R.matlab)

setwd("~/Dropbox/R/NPFDA")

library(termstrc)
library(tseries)
library(forecast)
source("biblioteca-procedimento-npfda.r")
intervalo <- c(1,256)
dir <- "./Resultados/22-jan/"


# intervalo.treino <- 1:2500
# intervalo.teste <- 2501:5000
# intervalo.validacao <- 5001:6530

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
# 3  6  9	12	15	18	21	24	30	36	48	60	72	84	96	108	120	144	180	240
# maturidades <- c(0.25 , 0.5,  0.75,  1,	1.25,	1.5,	1.75,	2,	2.5,	3,	4,	5,	6,	7,	8,	9,	10,	12,	15, 20)
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,	1.5,	1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10,	12, 15, 20)
percentual.testar <- 0.2
intervalo.base <- 1:6535
taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [intervalo.base ,2:(length(tempo.maturidades)+1)])
datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [intervalo.base, 1]), format = "%Y-%m-%d")

#### Inicio estimações
simulacoes <- vector("list", 0)
parametro <- 0.5
horizonte <- 125
betas <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro,
                               datas, parametro)


for (j in 1:length(intervalos)) {
  intervalo.total <- intervalos[[j]]
  intervalo.passado <- IntervaloPassado(intervalo.total,percentual.testar)
  intervalo.futuro <- IntervaloFuturo(intervalo.total,percentual.testar, horizonte)
  ## Prevê Diebold-Li ##
#   betas.previstos <- DieboldLi.PreveBetas(betas, intervalo.passado,
#                                           intervalo.futuro)
#   taxas.previstas.diebold <- DieboldLi.BetasParaTaxas(
#     betas.previstos, tempo.maturidades, parametro)
  for (i in 1:length(vetor.maturidade)) {
    maturidade <- vetor.maturidade[i]
    ## Prevê RW sem drift ##
#     taxas.previstas.rw <- taxas.juro[py.range(intervalo.futuro) - horizonte, maturidade]
    ## Prevê Arima ##
    taxas.previstas.ar <- Arima.Previsao(taxas.juro[ , maturidade], 
                                      intervalo.passado,intervalo.futuro)
    ## Prevê Corte ##
#     retirar <- taxas.juro[, 6]
#     curvas <- PreparaCurvasCorte(taxas.juro,percentual.testar,intervalo.total,s=horizonte,maturidade=maturidade, retirar = retirar)
#     semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
#     taxas.previstas.corte.deriv <- predict(curvas,semimetricas)
#     
#     semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
#     taxas.previstas.corte.pca <- predict(curvas,semimetricas)
    
    ## Valores reais ##
#     valores.reais <- taxas.juro[py.range(intervalo.futuro), maturidade]
    
    ## Grava no arquivo ##
#     simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]] <- taxas.previstas.diebold[, maturidade]
#     simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]][["1"]] <- taxas.previstas.rw
    simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]][["1"]] <- taxas.previstas.ar
#     simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["1"]] <- taxas.previstas.corte.deriv
#     simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]] <- taxas.previstas.corte.pca
#     simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]][["0"]] <- valores.reais
  }
}

save(simulacoes,file=paste0(dir,"simulacoes"))


PlotarMetodo <- function(arquivo,intervalo,horizonte,metodo,parametro) {
  Plot.list(simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]])
  for (i in simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]]) 
}

PlotarExperimento <- function(arquivo,intervalo,maturidade,horizonte) {

  eqm <- array(0,dim=c(length(intervalos),5))
  dimnames(eqm) <- list(intervalos,c("diebold-li","rw","arima","corte.deriv","corte.pca"))
  
  for (i in 1:length(intervalos)) {
    intervalo <- intervalos[i]
    Plot.list(arquivo[[nome.vetor(intervalo)]][[as.character(horizonte)]][[as.character(maturidade)]])
    title(nome.vetor(intervalo))
    taxas.diebold <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]]
    taxas.rw <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]][["1"]]
    taxas.ar <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]][["1"]]
    taxas.corte.deriv <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["1"]]
    taxas.corte.pca <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]]
    valores.reais <- arquivo[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]][["0"]]    
    eqm[i,1] <- EQM(valores.reais,taxas.diebold); eqm[i,2] <- EQM(valores.reais,taxas.rw)
    eqm[i,3] <- EQM(valores.reais,taxas.ar);eqm[i,4] <- EQM(valores.reais, taxas.corte.deriv)
    eqm[i,5] <- EQM(valores.reais, taxas.corte.pca)    
  }  
    print(eqm)
}


PlotarExperimento <- function(arquivo,maturidade,horizonte,intervalos) {
  
  eqm <- array(0,dim=c(length(intervalos),5))
  dimnames(eqm) <- list(intervalos,c("diebold-li","rw","arima","corte.deriv","corte.pca"))
  
  for (i in 1:length(intervalos)) {
    intervalo.total <- intervalos[[i]]
    Plot.list(arquivo[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]])
    title(nome.vetor(intervalo.total))
    taxas.diebold <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]]
    taxas.rw <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]][["1"]]
    taxas.ar <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]][["1"]]
    taxas.corte.deriv <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["1"]]
    taxas.corte.pca <- simulacoes[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]]
    valores.reais <- arquivo[[nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]][["0"]]    
    eqm[i,1] <- EQM(valores.reais,taxas.diebold); eqm[i,2] <- EQM(valores.reais,taxas.rw)
    eqm[i,3] <- EQM(valores.reais,taxas.ar);eqm[i,4] <- EQM(valores.reais, taxas.corte.deriv)
    eqm[i,5] <- EQM(valores.reais, taxas.corte.pca)    
  }  
  print(eqm)
}



PlotarExperimento(simulacoes,2,horizonte,intervalos)
Multiplot.ts(matrix(rep(1:20,each=20),nrow=20))
