# Funções --------------------------------------------------------------


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

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
# 3  6  9	12	15	18	21	24	30	36	48	60	72	84	96	108	120	144	180	240
# maturidades <- c(0.25 , 0.5,  0.75,  1,	1.25,	1.5,	1.75,	2,	2.5,	3,	4,	5,	6,	7,	8,	9,	10,	12,	15, 20)
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,	1.5,	1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10,	12, 15, 20)
percentual.testar <- 0.2

## A base de dados pode ter todas as observações ou ser de 5 em 5 dias
## 
taxas.juro.inteira <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [, 2:(length(tempo.maturidades)+1)])
taxas.juro.reduzida <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [1:1307*5, 2:(length(tempo.maturidades)+1)])
taxas.juro <- taxas.juro.reduzida

## 
datas.inteira <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [ , 1]), format = "%Y-%m-%d")
datas.reduzida <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [ , 1]), format = "%Y-%m-%d")
datas <- datas.reduzida
#intervalos <- list(c(1,4000),c(1001,5000)) # para base completa
intervalos <- list(c(1,800),c(201,1000)) # para base reduzida
vetor.maturidade = c(1:4,6,8,12,17)


#### Inicio estimações
simulacoes <- vector("list", 0) # reseta os dados de simulacoes
parametro <- 0.5 # lambda para diebold-li
horizonte <- 5
betas.02 <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro,
                                 datas, 0.2)
betas.05 <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro,
                                     datas, 0.5)

for (j in 1:length(intervalos)) {
  intervalo.total <- intervalos[[j]]
  intervalo.passado <- IntervaloPassado(intervalo.total,percentual.testar)
  intervalo.futuro <- IntervaloFuturo(intervalo.total,percentual.testar, horizonte)
  ## Prevê Diebold-Li ##
  betas.previstos.05 <- DieboldLi.PreveBetas(betas.05, intervalo.passado,
                                          intervalo.futuro)
  taxas.previstas.diebold.05 <- DieboldLi.BetasParaTaxas(
    betas.previstos.05, tempo.maturidades, 0.5)
  betas.previstos.02 <- DieboldLi.PreveBetas(betas.02, intervalo.passado,
                                          intervalo.futuro)
  taxas.previstas.diebold.02 <- DieboldLi.BetasParaTaxas(
    betas.previstos.02, tempo.maturidades, 0.2)  
  for (i in 1:length(vetor.maturidade)) {
    maturidade <- vetor.maturidade[i]
    ## Prevê RW sem drift ##
#     taxas.previstas.rw <- taxas.juro[py.range(intervalo.futuro) - horizonte, maturidade]
    ## Prevê Arima ##
#     taxas.previstas.ar <- Ar.Previsao(taxas.juro[ , maturidade], 
#                                       intervalo.passado,intervalo.futuro)
    ## Prevê Corte ##
#     retirar <- taxas.juro[, 6]
#     curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro,retirar = retirar)
#     semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
#     taxas.previstas.corte.deriv <- predict(curvas,semimetricas)
#     semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
#     taxas.previstas.corte.pca <- predict(curvas,semimetricas)
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["1"]] <- taxas.previstas.corte.deriv
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]] <- taxas.previstas.corte.pca
#     
#     retirar <- betas.05[, 1]
#     curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro,retirar = retirar)
#     semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
#     taxas.previstas.corte.deriv <- predict(curvas,semimetricas)    
#     semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
#     taxas.previstas.corte.pca <- predict(curvas,semimetricas) ; rm(semimetricas)
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["beta1"]] <- taxas.previstas.corte.deriv
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["beta5"]] <- taxas.previstas.corte.pca
    ## Valores reais ##
#     valores.reais <- taxas.juro[py.range(intervalo.futuro), maturidade]
    
    ## Grava no arquivo ##
    simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]] <- taxas.previstas.diebold.05[, maturidade]
    simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.2"]] <- taxas.previstas.diebold.02[, maturidade]
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]] <- taxas.previstas.rw
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]] <- taxas.previstas.ar
#     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]] <- valores.reais
  }
}

# save(simulacoes,file=paste0(dir,"simulacoes"))
save(simulacoes,file="simulacoes_reduzida.dad")

PlotarMetodo <- function(arquivo,intervalo,horizonte,metodos,maturidade,parametro = NULL) {
###########################################
# Plota os métodos desejados
###########################################
  for(i in 1:length(metodos))
    
  if (missing(parametro)) {
    Plot.list(arquivo[[Nome.vetor(intervalo)]][[as.character(horizonte)]][[as.character(maturidade)]][[as.character(metodo)]])
  }
}

PlotarExperimento <- function(arquivo,maturidade,horizonte,intervalos) {
######################################################################
# Plota a série temporal para cada um dos métodos diferentes. A ordem está na
# variável "metodos". As cores para cada método estão na mesma ordem.
######################################################################
  metodos <- c("diebold-li","rw","arima","corte.deriv","corte.pca","combinado")
  eqm <- array(0,dim=c(length(intervalos),length(metodos)))
  dimnames(eqm) <- list(intervalos,metodos)
  
  for (i in 1:length(intervalos)) {
    intervalo.total <- intervalos[[i]]
    Plot.list(arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]])
    title(paste(Nome.vetor(intervalo.total)," | hor ", as.character(horizonte), " | maturidade", as.character(maturidade)))
    taxas.diebold <- simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]]
    taxas.rw <- simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]][["1"]]
    taxas.ar <- simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]][["1"]]
    taxas.corte.deriv <- simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["1"]]
    taxas.corte.pca <- simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]]
    taxas.combinado <- (taxas.corte.pca + taxas.corte.deriv)/2
    lines.ts(taxas.combinado, col = 7)
    valores.reais <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]][["0"]]
    eqm[i,1] <- EQM(valores.reais,taxas.diebold); eqm[i,2] <- EQM(valores.reais,taxas.rw)
    eqm[i,3] <- EQM(valores.reais,taxas.ar);eqm[i,4] <- EQM(valores.reais, taxas.corte.deriv)
    eqm[i,5] <- EQM(valores.reais, taxas.corte.pca);eqm[i,6] <- EQM(valores.reais, taxas.combinado)
  }  
  print(eqm)
}

PlotarExperimentoErros <- function(arquivo,maturidade,horizonte,intervalos) {
######################################################################
# Plota a série temporal de erros para cada um dos métodos diferentes. A ordem está na
# variável "metodos". As cores para cada método estão na mesma ordem.
######################################################################  
  metodos <- c("diebold-li","rw","arima","corte.deriv","corte.pca","combinado","corte.deriv.beta","corte.pca.beta")
  eqm <- array(0,dim=c(length(intervalos),length(metodos)))
  dimnames(eqm) <- list(intervalos,metodos)
  for (i in 1:length(intervalos)) {
    intervalo.total <- intervalos[[i]]
    taxas.diebold <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]]
    taxas.rw <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]][["1"]]
    taxas.ar <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]][["1"]]
    taxas.corte.deriv <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["1"]]
    taxas.corte.pca <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]]
    taxas.combinado <- (taxas.corte.pca + taxas.corte.deriv)/2
    taxas.corte.deriv <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv"]][["beta1"]]
    taxas.corte.pca <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["beta5"]]
    #lines.ts((taxas.combinado-series.erros[, dim(series.erros)[2]])^2, col = 7)
    valores.reais <- arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]][["0"]]    
    eqm[i,1] <- EQM(valores.reais,taxas.diebold); eqm[i,2] <- EQM(valores.reais,taxas.rw)
    eqm[i,3] <- EQM(valores.reais,taxas.ar);eqm[i,4] <- EQM(valores.reais, taxas.corte.deriv)
    eqm[i,5] <- EQM(valores.reais, taxas.corte.pca);eqm[i,6] <- EQM(valores.reais, taxas.combinado)
    series <- as.data.frame(arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]])
    series <- cbind(series,taxas.combinado)
    series.erros <- (series-valores.reais)^2
    Multiplot.ts(series.erros)
    # Plot.list(arquivo[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]])
    title(paste("Série erros: ",Nome.vetor(intervalo.total)," | hor ", as.character(horizonte), " | maturidade", as.character(maturidade)))
  }  
  print(eqm)
}

Tabelao <- function(arquivo,vetor.maturidades,vetor.horizontes,vetor.intervalos) {
#######################################################
#
#
######################################################
for (i in vetor.maturidades)
  for (j in vetor.horizontes)
    for(k in vetor.intervalos) {
      series <- as.data.frame(arquivo[[Nome.vetor(k)]][[as.character(j)]][[as.character(i)]])
      eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
      print(paste(i," ",j," ",k[1],"-",k[2],":",sep='')); print(eqm)
    }
}

Combinar.Previsoes <- function(df,m1,m2) {
####################
  df <- cbind(df,(df[, which(names(series)==m1)] + df[, which(names(series)==m2)])/2)
  names(df)[ncol(df)] <- paste(m1,"+",m2)
  return(df)
}


Plotar2 <- function(arquivo,maturidade,horizonte,intervalo) {
  #######################################################
  # 
  #
  ######################################################
  series <- as.data.frame(arquivo[[Nome.vetor(intervalo)]][[as.character(horizonte)]][[as.character(maturidade)]])
  series <- Combinar.Previsoes(series,"corte.pca.5","diebold.0.5")
  series.erros <- (series[, -which(names(series)=="valores.reais")] - series[, which(names(series)=="valores.reais")])^2
  series.erros.medios <- apply(series.erros,2,FUN=cummean)
  eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
  Multiplot.ts(series); title("Series")
  Multiplot.ts(series.erros); title("Series dos erros")
  Multiplot.ts(series.erros.medios); title("Series dos erros médios")
  print(eqm)
}


# DESCOBRIR A ORDEM DAS CORES
# plot(1:10,rep(0,10),col=1:10,lwd=5)

curvas <- PreparaCurvasCorte(taxas.juro,percentual.testar,intervalo=c(1,2500),s=125,maturidade=17,retirar=betas[, 1])
semi <- SemimetricasClasse(curvas)
previsto <- predict(curvas,semi)
IntervaloFuturo(c(1,2500),percentual.testar, horizonte)
valores.reais <- taxas.juro[py.range(intervalo.futuro), 17]

EQM(valores.reais,previsto)
lines.ts(previsto,col=8)
Plotar2(simulacoes,maturidade=2,horizonte=5,c(1,800))
PlotarExperimentoErros(simulacoes,17,horizonte,list(c(1,800)))
Multiplot.ts(matrix(rep(1:20,each=20),nrow=20))
View(as.data.frame(simulacoes[[Nome.vetor(intervalo.total)]][[as.character(25)]][[as.character(maturidade)]]))
Tabelao(simulacoes,vetor.maturidade,c(5,25,50),list(c(1,800),c(201,1000)))
