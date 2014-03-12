# Código principal --------------------------------------------------------

load(".RData")
setwd("~/Dropbox/R/NPFDA")

library(termstrc)
source("Experimentos reproduziveis/biblioteca-diebold-li.R")
source("Experimentos reproduziveis/biblioteca-npfda-ettj.R")
source("Experimentos reproduziveis/biblioteca simulacoes.R")
source("biblioteca-npfda.r")

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,  1.5,	1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10,	12, 15, 20) ## Usar para Base nefasta
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,  1.5,  1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10) ## Usar para Base mensal


# Quanto do intervalo total será usado para testar; 0.5 = 50% É utilizado para 
# definir as variáveis intervalo.passado e intervalo.futuro. Essas duas 
# variáveis são necessárias para dar entrada nas funções de previsão. 
# intervalo.passado é um vetor de dimensão 2 que indica o início e o fim do 
# intervalo que será utilizado como treino. Já intervalo.futuro também é um 
# vetor de dimensão 2, indicando o início e o fim do intervalo que a previsão 
# será feita. Exemplo: sendo o valor das duas variáveis, respectivamente, 
# c(1,1000) e c(1200,1400), a estimação é realizada no intervalo 1:1000 e 201 
# valores serão previstos, para o intervalo 1200:1400. Finalmente, o que a 
# variável percentual testar faz, portanto, é definir o tamanho relativo entre o
# intervalo.passado e intervalo.futuro. Consulte também a documentação das
# funções IntervaloPassado e IntervaloFuturo.
percentual.testar <- 0.2



############## Escolher UMA das BASES DE DADO abaixo: ######################
## 1.a) BASE NEFASTA: A base de dados pode ter todas as observações ...
taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [, 2:(length(tempo.maturidades)+1)])
datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [ , 1]), format = "%Y-%m-%d")
## 1.b)  ... ou conter as observações de 5 em 5 dias
taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [1:1307*5, 2:(length(tempo.maturidades)+1)])
datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [1:1307*5 , 1]), format = "%Y-%m-%d")
## 2) Base menor
taxas.juro <- as.matrix(read.table("~/Dropbox/R/NPFDA/Dados/dados_mensais.csv", 
                                   header=T, quote="\""))
datas <- as.Date(as.character(rownames(taxas.juro)), format = "%Y%m%d")
##################################

# intervalos: uma lista com os intervalos nos quais se fará a análise
vetor.intervalos <- list(c(1,800),c(101,900),c(201,1000),c(301,1100),
                         c(401,1200),c(501,1300)) # para base reduzida

# vetor.maturidade: contém os índices das maturidades que serão estimadas 
# e previstas
vetor.maturidade = c(1:4,6,8,12,17)

# vetor.horizontes
vetor.horizontes <- c(5,25,50)  vetor.horizontes <- c(1,3,6,12)

#### Inicio estimações
# simulacoes <- vector("list", 0) # RESETA os dados de simulacoes
# load("simulacoes_reduzida.dad")
# load("simulacoes_reduzida_base_mensal.dad")

betas.02.a <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro + 3,
                                  datas, 0.2)
betas.05.a <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro + 3,
                                  datas, 0.5)
betas.02 <- cbind(beta_1=betas.02.a[, 1] - 3,betas.02.a[, 2:3]); rm(betas.02.a)
betas.05 <- cbind(beta_1=betas.05.a[, 1] - 3,betas.05.a[, 2:3]); rm(betas.05.a)

for (horizonte in vetor.horizontes) {
  for (j in 1:length(vetor.intervalos)) {
    intervalo.total <- vetor.intervalos[[j]]
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
      ## Grava previsões de diebold-li para as maturidades  desejadas ##
            simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold.05"]] <- taxas.previstas.diebold.05[, maturidade]
            simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold.02"]] <- taxas.previstas.diebold.02[, maturidade]
      ## Prevê RW sem drift ##
          taxas.previstas.rw <- taxas.juro[py.range(intervalo.futuro) - horizonte, maturidade]
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]] <- taxas.previstas.rw
      ## Prevê Arima ##
          taxas.previstas.ar <- Ar.Previsao(taxas.juro[ , maturidade], 
                                            intervalo.passado,intervalo.futuro)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]] <- taxas.previstas.ar
      ## Prevê Corte ##
          retirar <- taxas.juro[, 6]
          curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro,retirar = retirar)
          semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
          taxas.previstas.corte.deriv.1 <- predict(curvas,semimetricas)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.d1.rmat6"]] <- taxas.previstas.corte.deriv.1
          semimetricas <- SemimetricasClasse(curvas,q=0,tipo="deriv")
          taxas.previstas.corte.deriv.0 <- predict(curvas,semimetricas)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.d0.rmat6"]] <- taxas.previstas.corte.deriv.0
          semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
          taxas.previstas.corte.pca <- predict(curvas,semimetricas)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.p5.rmat6"]] <- taxas.previstas.corte.pca

          curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro)
          semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
          taxas.previstas.corte.deriv <- predict(curvas,semimetricas)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.d1"]] <- taxas.previstas.corte.deriv
          semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
          taxas.previstas.corte.pca <- predict(curvas,semimetricas)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.p5"]] <- taxas.previstas.corte.pca
          
          retirar <- betas.05[, 1]
          curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro,retirar = retirar)
          semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
          taxas.previstas.corte.deriv <- predict(curvas,semimetricas)    
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.d1.rb1"]] <- taxas.previstas.corte.deriv
          semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
          taxas.previstas.corte.pca <- predict(curvas,semimetricas)
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["crt.p5.rb1"]] <- taxas.previstas.corte.pca
          rm(semimetricas)
      ## Valores reais ##
          valores.reais <- taxas.juro[py.range(intervalo.futuro), maturidade]
          simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]] <- valores.reais      
    }
  }
}

# salva o objeto simulações no lugar desejado. Aconselha-se salvar com
# frequência para não correr o risco de sobrescrever os experimentos e perder
# todo o trabalho
# save(simulacoes,file="simulacoes_reduzida.dad")
save(simulacoes,file="simulacoes_reduzida_base_mensal.dad")

# DESCOBRIR A ORDEM DAS CORES
# plot(1:10,rep(0,10),col=1:10,lwd=5)

curvas <- PreparaCurvasCorte(taxas.juro,percentual.testar,intervalo=c(1,2500),s=125,maturidade=17,retirar=betas[, 1])
semi <- SemimetricasClasse(curvas)
previsto <- predict(curvas,semi)
IntervaloFuturo(c(1,2500),percentual.testar, horizonte)
valores.reais <- taxas.juro[py.range(intervalo.futuro), 17]

EQM(valores.reais,previsto)
lines.ts(previsto,col=8)
PlotarSeries(simulacoes,maturidade=2,horizonte=5,intervalo=c(201,1000))
PlotarExperimentoErros(simulacoes,17,horizonte,list(c(1,800)))
Multiplot.ts(matrix(rep(1:20,each=20),nrow=20))
View(as.data.frame(simulacoes[[Nome.vetor(intervalo.total)]][[as.character(25)]][[as.character(maturidade)]]))


tabelao <- Tabelao(simulacoes,2,c(5,25,50),vetor.intervalos)
tabelao <- tabelao[, c(-8,-9)] # tira as estimacoes do que não retiram nada
require(stringr)
par.ordenar <- 1; tabelao <- tabelao[order(str_split_fixed(rownames(tabelao), "[|]",3)[, par.ordenar]), ] # par.ordenar define qual por qual dos 3 campos que será ordenada as estimações
# tabelao <- CombinarPrevisoes(tabelao,"corte.pca.5","rw")
# tabelao <- CombinarPrevisoes(tabelao,"corte.deriv.1","rw")
heatmap(tabelao, scale="row",Rowv=NA) # scale="row",Rowv=NA
require(lattice)
levelplot(tabelao, scale=list(x=list(rot=45)))
summary(tabelao)

tabelao <- Tabelao(simulacoes,1,5,vetor.intervalos)
heatmap(tabelao[, c(-8,-9)],Rowv=NA) # tira as estimacoes do que não retiram nada
levelplot(t(tabelao[order(rownames(tabelao)), c(-8,-9)]), scale=list(x=list(rot=45)))
View(tabelao)
tabelao.agrupado <- TabelaoAgrupado(simulacoes,vetor.maturidade,horizonte=50,vetor.intervalos)
heatmap(tabelao.agrupado[c(-8,-9), ],scale="column",Colv=NA) # tira as estimacoes do que não retiram nada

PlotarSeries(simulacoes,maturidade=4,horizonte=25,intervalo=c(1,800))
View(tabelao.agrupado)

series <- ExtraiSeries(simulacoes,1,5,c(1,800))
gw.test(x = tmp[,7], y = tmp[,3], p = tmp[,12], T = 800, tau = 5,method="HAC")
GiacominiWhite(series,benchmarking=3,valor.real=12,tamanho=800,horizonte=5)
View(series)
