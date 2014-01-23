
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


vetor.truncar <- list(c(1,1000),c(201,1200),c(401,1400),c(601,1600)) # ,c(201,1200),c(401,1400),c(601,1600)
vetor.grau.derivada = c(0,1,2) #quais derivadas testar 
vetor.n.componentes.principais = c(5,6,7) # quais componentes principais testar
vetor.s = c(21,125)
vetor.obs.por.curva = c( 50) # 30,40, 50, 60,
vetor.n.vizinhos = c(1:10*2)
vetor.maturidade = c(1:4,6,8,12,17)
vetor.passo = c(5,15,20,30)


# Código principal --------------------------------------------------------

#library(R.matlab)

setwd("~/Dropbox/R/NPFDA")

source("biblioteca-procedimento-npfda.r")
#base <- read.csv("~/Dropbox/R/NPFDA/base_nefasta_mat.csv", sep=",")[2:21]
base <- read.csv("Dados/base_nefasta_mat.csv", sep=",")[2:21]

temp <- read.table("~/Dropbox/R/NPFDA/base_nefasta_mat.csv", sep=",",header=TRUE)
z <- zoo(temp[, 2:21], as.Date(as.character(temp[, 1]), format = "%Y-%m-%d"))

base.corte <- base - base[, 1]
plot(maturidades,base[1000, ],type="o")

EstimacaoCompletaPasso(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s,vetor.passo,vetor.obs.por.curva,vetor.n.componentes.principais,vetor.grau.derivada)
EstimacaoCompletaRW(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s)
EstimacaoCompletaCorte(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s,todas.maturidades,vetor.n.componentes.principais,vetor.grau.derivada)

Multiplot.ts(base)
