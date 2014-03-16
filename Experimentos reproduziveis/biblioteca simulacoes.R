###########################################################
# Conjunto de funções 
#
#
###########################################################


AnimacaoJuros <- function(taxas.juro, tempo.maturidades, file) {
  
  image <- array(0, c(2,17, dim(taxas.juro)[1]))
  for (i in 1:(dim(taxas.juro)[1])) {
    image[,,i] <- as.matrix(cbind(tempo.maturidades,taxas.juro[i,]))
  }
  View(image[,,i])
}

Ar.Previsao <- function(ts, intervalo.passado, intervalo.futuro) {
  #########################################################
  # Realiza previsões para o intervalo futuro. O horizonte é definido pela
  # diferença entre os intervalos e o modelo é um AR(1) por default
  # ARGS
  #   ts: a série temporal para prever
  #   intervalo.passado: vetor de tamanho 2 com o início e fim do intervalo de treino
  #   intervalo.futuro: vetor de tamanho 2 com o início e fim do intervalo para previsão
  #   grau.arima: vetor de ordem 3 que especifica os graus do modelo arima c(p,i,q)
  # RETORNA
  # Um vetor com as previsões
  #########################################################
  previstos <- rep(0,times = diff(intervalo.futuro) + 1)
  horizonte  <- intervalo.futuro[1] - intervalo.passado[2]
  for (i in intervalo.futuro[1]:intervalo.futuro[2]) {
    estimacao.arima <- ar(ts[intervalo.passado[1]:(i - horizonte)],order.max=1,method="ols")
    previstos[i - intervalo.futuro[1] + 1] <- predict(estimacao.arima,
                                                      n.ahead=horizonte, aic = FALSE)$pred[horizonte]
  }
  return(previstos)
}



CombinarPrevisoes <- function(df, ...) {
  ####################
  df <- as.data.frame(df)
  metodos <- c(...)
  tmp <- apply(sapply(metodos, function (x) df[, which(names(df)==x)]),1,sum)/
    length(metodos)
  df <- cbind(tmp, df)
  names(df)[1] <- paste(metodos,collapse=' + ')
  return(as.matrix(df))
}


cummean <- function(x) {
  ############################################################ 
  #Calculates the cumulative mean of a vector. Works like the function cumsum, but
  #taking its mean
  ###########################################################
  if (is.vector(x)) output <- cumsum(x)/1:length(x)
  return(output)
}

EQM <- function(d1,d2) {
  ##########################################################
  # Calcula o erro quadrático médio entre dois vetores
  #########################################################
  mean((d1-d2)^2)
}


ExtraiSeries <- function(arquivo,maturidade,horizonte,vetor.intervalos) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
  # de maturidade x horizonte x intervalo.
  ######################################################
  series <- NULL
  for (intervalo in vetor.intervalos) {
    series.tmp <- as.data.frame(arquivo[[Nome.vetor(intervalo)]][[as.character(horizonte)]][[as.character(maturidade)]])
    series <- rbind(series,series.tmp)
  }
  return(series)
}


GiacominiWhite <- function(arquivo, benchmarking, valor.real, tamanho, horizonte, vetor.maturidade,vetor.intervalos) {
  ############################################################
  # Aplica o teste de Giacomini-White para uma sequência de previsões, comparando
  # todas com a capacidade de previsão de um mesmo método como benchmarking
  # ARGS
  #   tabela: tabela com os valores previstos para cada um dos métodos, incluindo
  #                 o benchmarking e os valores.reais. Séries nas colunas.
  #   benchmarking: o índice da coluna da variável para tomar como benchmarking
  #   valor.real: o índice da coluna dos valores da série real
  #   tamanho: o tamanho da amostra que gerou as previsões
  #   horizonte: o horizonte das previsões
  #   vetor.maturidade: vetor com o índice das maturidades para testar
  # RETORNO
  #   um vetor com os valores da estatística de Giacomini-White, para cada  método
  ############################################################# 
  #install.packages(c("polynom", "fracdiff", "hypergeo", "longmemo")) 
  #install.packages("~/Dropbox/R/NPFDA/afmtools_0.1.8.tar.gz", repos = NULL, type
  #= "source")
  tabela.gw <- NULL
  for (i in vetor.maturidade) {
    tabela <- ExtraiSeries(arquivo,maturidade=i,horizonte=5,vetor.intervalos)
    if (!require(afmtools)) stop("Por favor, instale a biblioteca afmtools. Instrucoes na descricao da funcao.")
    tmp <- apply(tabela[, c(-benchmarking,-valor.real)], 2, gw.test, 
                 y = tabela[, benchmarking], p = tabela[, valor.real],
                 T = tamanho, tau = horizonte, method="HAC")
    resultados <- unlist(lapply(tmp,function(x) x$p.value))
    #gw.test(x = tmp[,7], y = tmp[,3], p = tmp[,12], T = 800, tau = 5,method="HAC")
    tabela.gw <- rbind(tabela.gw,i=resultados)
  }
  return(tabela.gw)
}


IntervaloFuturo <- function(intervalo,percentual.testar, horizonte) {
  #############################################################
  # ARGS
  #   intervalo: vetor tipo c(x,y)
  #   percentual.testar: um número entre 0 e 1
  # RETORNO
  #   um vetor c(z,y), em que z é o primeiro valor que será previsto  
  #   Ex: IntervaloFuturo(c(1001,2000),0.8,h) retorna c(1801 + h,2000)
  #############################################################
  
  tamanho.serie.temporal  = intervalo[2]-intervalo[1]+1
  output <- c(trunc((1-percentual.testar)*tamanho.serie.temporal
  ) + intervalo[1]-1 + horizonte, tamanho.serie.temporal + 
    intervalo[1] - 1)
  return(output)
}

IntervaloPassado <- function(intervalo,percentual.testar) {
  #############################################################
  # ARGS
  #   intervalo: vetor tipo c(x,y)
  #   percentual.testar: um número entre 0 e 1
  # RETORNO
  #   um vetor c(x,z), em que z corresponde ao (1-percentual.testar)
  #   do intervalo inicial. 
  #   Ex: Intervalo.base.previsao(c(1001,2000),0.8) retorna c(1001,1800)
  #############################################################
  output <- intervalo
  output[2] <- trunc((intervalo[2]-intervalo[1]+1)*
                       (1-percentual.testar) + intervalo[1] - 1)
  return(output)
}




Multiplot.ts <- function(data, ... , range, dates){
  ############################################################
  # This function plots many time series in the same graph
  #
  # ARGS  
  #    Data is a matrix with each serie in its lines
  #    Range is an optional vector with the indexes of the time series
  #       to be exibited. If it is absent, then every line is ploted
  ############################################################
  if (missing(range)) 
    range <- 1:dim(data)[2]
  size.series <- dim(data)[2]
  if (missing(dates)) {
    plot(data[, range[1]],col=1,ylim=c(min(data[, range]), 
                                       max(data[, range])), type="l",  ...)    
    for (i in range[-1]) {
      lines(data[, i],col=i, lty=((i - 1) %/% 8 + 1), ...)   
    }
  } else {
    plot(dates, data[, range[1]],col=1,ylim=c(min(data[, range]), 
                                              max(data[, range])), type="l",  ...)
    for (i in range[-1]) {
      lines(dates,data[, i],col=i, lty=((i - 1) %/% 8 + 1), ...)   
    }
  } 
}  

Nome.vetor <- function(object) {
  return(paste0(object[1],":",object[2]))
}

Plot.list <- function(object) {
  tmp <- as.data.frame(object)
  Multiplot.ts(tmp)
}

PlotarSeries <- function(arquivo,percentual.testar,maturidade,horizonte,subconjunto = NULL) {
  ####################################################### 
  # Exibe 3 gráficos: I) séries previstas por cada método e a realizada; II) erro
  # quadrático de cada série prevista, por unidade de tempo; III) erro quadrático
  # médio de cada série, até a unidade de tempo t. 
  ######################################################
  series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
  if (!is.null(subconjunto)) series <- series[py.range(subconjunto), ]
  series <- series[, order(names(series))] # coloca as colunas em ordem alfabética
  #series <- Combinar.Previsoes(series,"corte.pca.5","diebold.0.5") # para combinar previsões
  series.erros <- as.data.frame((series[, -which(names(series)=="valores.reais")] - series[, which(names(series)=="valores.reais")])^2)
  series.erros.medios <- as.data.frame(apply(series.erros,2,FUN=cummean))
  eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
  parte.titulo <- paste("horizonte ", as.character(horizonte), " | maturidade", as.character(maturidade))
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

Tabelao <- function(arquivo,vetor.maturidades,vetor.horizontes,vetor.intervalos) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
  # de maturidade x horizonte x intervalo.
  ######################################################
  tabela <- NULL
  for (i in vetor.maturidades)
    for (j in vetor.horizontes)
      for(k in vetor.intervalos) {
        series <- as.data.frame(arquivo[[Nome.vetor(k)]][[as.character(j)]][[as.character(i)]])
        eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
        nome.linha <- paste(j,"|",tempo.maturidades[i],"|",k[1],"-",k[2],sep='')
        tabela <- rbind(tabela,eqm)
        rownames(tabela)[dim(tabela)[1]] <- nome.linha
      }
  return(tabela)
}

TabelaoAgrupado <- function(arquivo,vetor.maturidade,horizonte,vetor.intervalos) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
  # de maturidade x horizonte x intervalo.
  ######################################################
  tabela <- NULL
  for (i in vetor.maturidade) {
    tabela.tmp <- Tabelao(arquivo, i, horizonte, vetor.intervalos)
    tabela <- cbind(tabela,colMeans(tabela.tmp))
  }
  colnames(tabela) <- tempo.maturidades[vetor.maturidade]
  return(t(tabela))
}

TabelaEQM <- function(arquivo,vetor.maturidade, horizonte, normaliza = TRUE, subconjunto = NULL) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada
  # de maturidade x método, para um determinado horizonte.
  ######################################################
  tabela <- NULL
  for (maturidade in vetor.maturidade) { 
    # arquivo <- simulacoes
    series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
    if (!is.null(subconjunto)) series <- series[py.range(subconjunto), ]
    eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
    tabela <- cbind(tabela,eqm)
  }
  colnames(tabela) <- tempo.maturidades[vetor.maturidade]
  if (normaliza) tabela <- tabela[-which(rownames(tabela)=="rw"), ] / tabela[which(rownames(tabela)=="rw"), ]
  return(t(tabela))
}
