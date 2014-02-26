# Código principal --------------------------------------------------------

setwd("~/Dropbox/R/NPFDA")

library(termstrc)
source("Experimentos reproduziveis/biblioteca-diebold-li.R")
source("Experimentos reproduziveis/biblioteca-npfda-ettj.R")
source("Experimentos reproduziveis/biblioteca simulacoes.R")
source("biblioteca-npfda.r")

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,  1.5,	1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10,	12, 15, 20)

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
## 1) A base de dados pode ter todas as observações ...
taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [, 2:(length(tempo.maturidades)+1)])
datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [ , 1]), format = "%Y-%m-%d")
## 2)  ... ou conter as observações de 5 em 5 dias
taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [1:1307*5, 2:(length(tempo.maturidades)+1)])
datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [1:1307*5 , 1]), format = "%Y-%m-%d")
##################################

# intervalos: uma lista com os intervalos nos quais se fará a análise
vetor.intervalos <- list(c(1,800),c(201,1000)) # para base reduzida

# vetor.maturidade: contém os índices das maturidades que serão estimadas 
# e previstas
vetor.maturidade = c(1:4,6,8,12,17)

# vetor.horizontes
vetor.horizontes <- c(5,25,50)

#### Inicio estimações
# simulacoes <- vector("list", 0) # RESETA os dados de simulacoes
load("simulacoes_reduzida.dad")


betas.02 <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro,
                                  datas, 0.2)
betas.05 <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro,
                                  datas, 0.5)
for (horizonte in vetor.horizontes) {
  for (j in 1:length(vetor.intervalos)) {
    intervalo.total <- vetor.intervalos[[j]]
    intervalo.passado <- IntervaloPassado(intervalo.total,percentual.testar)
    intervalo.futuro <- IntervaloFuturo(intervalo.total,percentual.testar, horizonte)
    ## Prevê Diebold-Li ##
#     betas.previstos.05 <- DieboldLi.PreveBetas(betas.05, intervalo.passado,
#                                                intervalo.futuro)
#     taxas.previstas.diebold.05 <- DieboldLi.BetasParaTaxas(
#       betas.previstos.05, tempo.maturidades, 0.5)
#     betas.previstos.02 <- DieboldLi.PreveBetas(betas.02, intervalo.passado,
#                                                intervalo.futuro)
#     taxas.previstas.diebold.02 <- DieboldLi.BetasParaTaxas(
#       betas.previstos.02, tempo.maturidades, 0.2)  
    for (i in 1:length(vetor.maturidade)) {
      maturidade <- vetor.maturidade[i]
      ## Prevê RW sem drift ##
      #     taxas.previstas.rw <- taxas.juro[py.range(intervalo.futuro) - horizonte, maturidade]
      ## Prevê Arima ##
      #     taxas.previstas.ar <- Ar.Previsao(taxas.juro[ , maturidade], 
      #                                       intervalo.passado,intervalo.futuro)
      ## Prevê Corte ##
            retirar <- taxas.juro[, 6]
            curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro,retirar = retirar)
#             semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
#             taxas.previstas.corte.deriv.1 <- predict(curvas,semimetricas)
#             simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv.0"]] <- taxas.previstas.corte.deriv.1
            semimetricas <- SemimetricasClasse(curvas,q=0,tipo="deriv")
            taxas.previstas.corte.deriv.0 <- predict(curvas,semimetricas)
            simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv.0"]] <- taxas.previstas.corte.deriv.0
      #     semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
#           taxas.previstas.corte.pca <- predict(curvas,semimetricas)
#           simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca"]][["5"]] <- taxas.previstas.corte.pca

#           curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro)
#           semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
#           taxas.previstas.corte.deriv <- predict(curvas,semimetricas)
#           semimetricas <- SemimetricasClasse(curvas,q=5,tipo="pca")
#           taxas.previstas.corte.pca <- predict(curvas,semimetricas)
#           simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.deriv1.retira0"]] <- taxas.previstas.corte.deriv
#           simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["corte.pca5.retira0"]] <- taxas.previstas.corte.pca     
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
#       simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.5"]] <- taxas.previstas.diebold.05[, maturidade]
#       simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["diebold"]][["0.2"]] <- taxas.previstas.diebold.02[, maturidade]
      #     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["rw"]] <- taxas.previstas.rw
      #     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["ar"]] <- taxas.previstas.ar
      #     simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores reais"]] <- valores.reais
    }
  }
}

# salva o objeto simulações no lugar desejado. Aconselha-se salvar com
# frequência para não correr o risco de sobrescrever os experimentos e perder
# todo o trabalho
save(simulacoes,file="simulacoes_reduzida.dad")


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


tabelao <- Tabelao(simulacoes,vetor.maturidade,c(5,25,50),list(c(1,800)))
heatmap(tabelao[, c(-9,-10)],scale="none",Rowv=NA) # tira as estimacoes do que não retiram nada
require(lattice)
levelplot(t(tabelao[order(rownames(tabelao), c(-9,-10)]), scale=list(x=list(rot=45)))
tabelao <- Tabelao(simulacoes,vetor.maturidade,c(5,25,50),list(c(201,1000)))
heatmap(tabelao[, c(-9,-10)],Rowv=NA) # tira as estimacoes do que não retiram nada
levelplot(t(tabelao[order(rownames(tabelao)), c(-9,-10)]), scale=list(x=list(rot=45)))
