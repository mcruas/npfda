# Código principal --------------------------------------------------------

setwd("~/Dropbox/R/NPFDA")

library(termstrc)
source("Experimentos reproduziveis/biblioteca-diebold-li.R")
source("Experimentos reproduziveis/biblioteca-npfda-ettj.R")
source("Experimentos reproduziveis/biblioteca simulacoes.R")

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,	1.5,	1.75,	2,	2.5,	3,	4,
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

simulacoes <- vector("list", 0) # RESETA os dados de simulacoes


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
}

# salva o objeto simulações no lugar desejado. Aconselha-se salvar com
# frequência para não correr o risco de sobrescrever os experimentos e perder
# todo o trabalho
save(simulacoes,file="simulacoes_reduzida.dad")


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
# eqm <- matrix(nrow=length(vetor.maturidades)*length(vetor.intervalos)*length(vetor.horizontes),ncol=)
  tabela <- NULL
  for (i in vetor.maturidades)
    for (j in vetor.horizontes)
      for(k in vetor.intervalos) {
        series <- as.data.frame(arquivo[[Nome.vetor(k)]][[as.character(j)]][[as.character(i)]])
        eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
        nome.linha <- paste(i," ",j," ",k[1],"-",k[2],sep='')
        tabela <- rbind(tabela,eqm)
        rownames(tabela)[dim(tabela)[1]] <- nome.linha
      }
 return(tabela)
}

Combinar.Previsoes <- function(df, ...) {
####################
  metodos <- c(...)
  tmp <- apply(sapply(metodos, function (x) df[, which(names(df)==x)]),1,sum)/
                                                      length(metodos)
  df <- cbind(tmp, df)
  names(df)[1] <- paste(metodos,collapse=' + ')
  return(df)
}

Plotar2 <- function(arquivo,percentual.testar,maturidade,horizonte,intervalo) {
  #######################################################
  # 
  #
  ######################################################
  series <- as.data.frame(arquivo[[Nome.vetor(intervalo)]][[as.character(horizonte)]][[as.character(maturidade)]])
  series <- series[, order(names(series))] # coloca as colunas em ordem alfabética
  #series <- Combinar.Previsoes(series,"corte.pca.5","diebold.0.5")
  series.erros <- as.data.frame((series[, -which(names(series)=="valores.reais")] - series[, which(names(series)=="valores.reais")])^2)
  series.erros.medios <- as.data.frame(apply(series.erros,2,FUN=cummean))
  eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
  # intervalo.futuro <- IntervaloFuturo(intervalo,percentual.testar, horizonte)
  parte.titulo <- paste(Nome.vetor(intervalo.total)," | hor ", as.character(horizonte), " | maturidade", as.character(maturidade))
  par(xpd=T, mar=par()$mar+c(0,0,0,4)) # Diminui o espaço do gráfico à direita para ter mais espaço para a legenda
  Multiplot.ts(series.erros.medios); title(paste("Series dos erros médios: ",parte.titulo),cex.main=0.8)
  legend("topright", inset=c(-0.32,0), legend=names(series.erros.medios), lty = (1:length(names(series.erros.medios)) - 1) %/% 8 + 1, col=1:length(names(series.erros.medios)), title="Métodos",cex=0.6)
  Multiplot.ts(series.erros); title(paste("Series dos erros: ",parte.titulo),cex.main=0.8)
  legend("topright", inset=c(-0.32,0), legend=names(series.erros), lty = (1:length(names(series)) - 1) %/% 8 + 1, col=1:length(names(series.erros)), title="Métodos",cex=0.6)
  Multiplot.ts(series,ylab="taxas de juro"); title(paste("Series: ",parte.titulo),cex.main=0.8)
  legend("topright", inset=c(-0.32,0), legend=names(series), lty = (1:length(names(series)) - 1) %/% 8 + 1, col=1:length(names(series)), title="Métodos",cex=0.6)
  par(mar=c(5, 4, 4, 2) + 0.1) # retorna 'par' para seu valor original
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
Plotar2(simulacoes,maturidade=2,horizonte=5,intervalo=c(1,800))
PlotarExperimentoErros(simulacoes,17,horizonte,list(c(1,800)))
Multiplot.ts(matrix(rep(1:20,each=20),nrow=20))
View(as.data.frame(simulacoes[[Nome.vetor(intervalo.total)]][[as.character(25)]][[as.character(maturidade)]]))
Tabelao(simulacoes,vetor.maturidade,c(5,25,50),list(c(1,800),c(201,1000)))
