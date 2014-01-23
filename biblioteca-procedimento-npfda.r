# Essa biblioteca possui funções para realizar previsões de séries temporais
# 
# https://chasquebox.ufrgs.br/data/public/3e5a40063049c87e15f39b8d2c9726f9.php?lang=pt-br


### bibliotecas necessárias para a boa execução das rotinas
source("biblioteca-npfda.r")
library("stringr")

DiretorioResultados  <- function() setwd("~/Dropbox/R/NPFDA/Resultados")
DiretorioPrincipal <-  function() setwd("~/Dropbox/R/NPFDA")

Ar.Previsao <- function(ts, intervalo.passado, intervalo.futuro, grau.arima = c(1,0,0)) {
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
  previstos <- rep(0,times=intervalo.futuro[2] - intervalo.futuro[1]+1)
  horizonte  <- intervalo.futuro[1] - intervalo.passado[2]
  for (i in intervalo.futuro[1]:intervalo.futuro[2]) {
    estimacao.arima <- arima(ts[intervalo.passado[1]:(i - horizonte)],order=grau.arima)
    previstos[i - intervalo.futuro[1] + 1] <- forecast(estimacao.arima,h=horizonte)$mean[horizonte]
  }
  return(previstos)
}


DieboldLi.EstimaBetas <- function (maturidades, taxas.juro, datas, lambda) {
#################################################################
# Esta função estima os valores de beta para uma dada maturidade, taxa de
# juros,datas e lambda.
# Precisa da biblioteca termstrc!
# ARGS
#   betas: matriz com a série temporal dos betas nas colunas
#   maturidades: vetor com a duração das maturidades
#   lambda: o valor de lambda utilizado na estimação
# RETORNA
#   uma matriz com as os coeficientes beta dos fatores em função do tempo
################################################################  
  
  # cria um objeto zeroyields e depois estima-se 
  datazeroyields1 <- zeroyields(maturidades, taxas.juro, datas)
  dl_res <- estim_nss(dataset=datazeroyields1, method="dl", lambda=lambda)
  output <- param(dl_res)[[1]]
  return(output)
  
}


DieboldLi.PreveBetas <- function(betas, intervalo.passado, intervalo.futuro) {
  ############################################################
  # Utiliza 
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   intervalo.passado: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo do treino
  #   intervalo.futuro: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo a ser previsto
  #   horizonte: o horizonte de previsão
  # RETORNA
  #   uma matriz com as taxas de juros em função do tempo e 
  #   das maturidades
  ###########################################################
  
  # estima previsões para cada beta
  betas.previstos <- matrix(data=0,nrow=(intervalo.futuro[2]-
                              intervalo.futuro[1]+1),ncol=3)
  horizonte  <- intervalo.futuro[1] - intervalo.passado[2]
  for (i in intervalo.futuro[1]:intervalo.futuro[2]) {
    
    # betas.arima é uma lista com o arima aplicado para cada um dos betas.
    # O beta de h passos à frente é obtido através de um AR(1)
    betas.arima <- apply(betas[intervalo.passado[1]:(i - horizonte), ],
                         2,auto.arima)
    betas.previstos[i - intervalo.futuro[1] + 1, ] <- 
                    sapply(lapply(betas.arima,forecast, h=horizonte),
                                  function(lst) lst$mean[horizonte])
  }
  return(betas.previstos)
}


DieboldLi.BetasRW <- function(betas, intervalo.passado, intervalo.futuro, horizonte) {
  ############################################################
  # Utiliza 
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   intervalo.passado: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo do treino
  #   intervalo.futuro: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo a ser previsto
  #   horizonte: o horizonte de previsão
  # RETORNA
  #   uma matriz com as taxas de juros em função do tempo e 
  #   das maturidades
  ###########################################################
  betas.previstos <- betas[py.range(intervalo.futuro) - horizonte, ]
  return(t(betas.previstos))
}


DieboldLi.BetasParaTaxas <- function(betas,maturidades,lambda) {
  ############################################################
  # Transforma a série temporal dos coeficientes dos fatores em
  # série temporal das taxas de juros para as maturidades desejadas
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   maturidades: vetor com a duração das maturidades
  #   lambda: o valor de lambda utilizado na estimação
  # RETORNA
  #   uma matriz com as taxas de juros em função do tempo e 
  #   das maturidades
  ###########################################################
  output <-  t(apply(betas,1,spr_dl,maturidades,lambda))
  colnames(output) <- maturidades * 12
  return(output)
}


EscreveResumo <- function (maturidade,curvas,s,valores.reais,truncar, ...) {
################################################################
# Essa função escreve um resumo das estimaçoes nos arquivo recebidos
# ARGS
#   maturidade: 
#
#
# RETORNA
#
################################################################

  mse  <- mean((curvas$estimacao - valores.reais)^2)
  if (!missing(...)){
    parametros <- (...)
  } else {
    parametros <- NULL
  }
  data=date()
  #linha=str_c(paste("mat",maturidade,sep=""),paste("hor",s,sep=""),paste(truncar[1],"-",truncar[2],sep=""),mse,curvas$tipo,parametros,data,sep=',')
  arquivo_para_escrever <- paste("Resultados/",paste(paste("mat",maturidade,sep=""),paste("hor",s,sep=""),paste(truncar[1],"-",truncar[2],sep=""),sep="_"),".csv",sep="")
  linha <- str_c(data,mse,curvas$tipo,parametros,sep=',')
  
  #escreve o resumo no arquivo
  cat(linha,file=arquivo_para_escrever,sep="\n",append=TRUE)
  
  #coloca num arquivo separado os dados estimados
  #cat(data,curvas$tipo,curvas$estimacao,file=nome_arquivos$registro,sep=",",append=TRUE)
  #cat("\n",file=nome_arquivos$registro,append=TRUE)
  
}



EstimacaoSerieTemporal <- function(serie.temporal,percentual.testar,s,passo,obs.por.curva,vetor.n.componentes.principais,vetor.grau.derivada) {
  
  curvas  <- NULL
  
  tamanho.serie <- length(serie.temporal)
  truncar <- c(1,tamanho.serie)
  valores.reais  <- PreparaValorFuturo(serie.temporal,percentual.testar,truncar,s,maturidade)      
  curvas <- PreparaCurvasPasso(base,percentual.testar,truncar,s,maturidade,passo,obs.por.curva)
  
  quantidade.estimacoes  <- length(vetor.n.componentes.principais)+length(vetor.grau.derivada)
  previsoes <- matrix(rnorm(800),nrow=quantidade.estimacoes,ncol=percentual.testar*tamanho.serie)
  
  # Estima utilizando a semimétrica PCA
  for (n.componentes.principais in vetor.n.componentes.principais) {
    
    semimetrica1 <- SemimetricPCA(curvas$passado.learn.media,curvas$passado.learn.media,n.componentes.principais)
    semimetrica2 <- SemimetricPCA(curvas$passado.learn.media,curvas$passado.test.media,n.componentes.principais)
    
    curvas$estimacao  <-  NULL
    curvas$estimacao  <- FunopareKnn(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,
                                       PRED=curvas$passado.test.media,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                       SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
        + curvas$media.test
    
  } # fim da PCA
  
  
  # Estima utilizando a semimétrica derivada
  for (grau.derivada in vetor.grau.derivada){
    
    semimetrica1 <- SemimetricDeriv(curvas$passado.learn.media,curvas$passado.learn.media,grau.derivada,nknot=5,range.grid=c(0,1))
    semimetrica2 <- SemimetricDeriv(curvas$passado.learn.media,curvas$passado.test.media,grau.derivada,nknot=5,range.grid=c(0,1))
    
    curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,PRED=curvas$passado.test.media,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values + curvas$media.test
    for (n.vizinhos in vetor.n.vizinhos)
    {
      curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,
                                      PRED=curvas$passado.test.media,kind.of.kernel="quadratic",
                                      SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
        + curvas$media.test
      #EscreveResumo(maturidade,curvas,s,valores.reais,truncar,n.vizinhos,paste(passo,",",obs.por.curva,",deriv regr ",grau.derivada,sep=""))
      
    }
  } # fim da deriv
  plot(serie.temporal,type="l")

}


EstimacaoCompletaCorte <- function(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s,todas.maturidades,vetor.n.componentes.principais,vetor.grau.derivada){
# a partir de uma base de dados, a estimação pelo método "corte" é executada para cada diferente parâmetro
# 
#
  
  curvas  <- NULL
  
  for (truncar in vetor.truncar)
    for (maturidade in vetor.maturidade)
      for (s in vetor.s){
        
        valores.reais  <- PreparaValorFuturo(base,percentual.testar,truncar,s,maturidade)      
        curvas <- PreparaCurvasCorte(base,percentual.testar,truncar,s,maturidade,todas.maturidades)
        
        semimetrica1 <- SemimetricDeriv(curvas$passado.learn,curvas$passado.learn,0,nknot=4,range.grid=c(0,1))
        semimetrica2 <- SemimetricDeriv(curvas$passado.learn,curvas$passado.test,0,nknot=4,range.grid=c(0,1))
        
        #        curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,PRED=curvas$passado.test,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values
        #        print(mean((valores.reais - curvas$estimacao)^2))
        
        for (n.vizinhos in vetor.n.vizinhos)
        {
          curvas$estimacao  <- FunopareKnn(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,
                                           PRED=curvas$passado.test,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                           SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values
          
          EscreveResumo(maturidade,curvas,s,valores.reais,truncar,paste("nviz ",n.vizinhos,"-deriv regr ",grau.derivada,sep=""))
          
        }
      } # fim da deriv

  semimetrica1 <- SemimetricPCA(curvas$passado.learn,curvas$passado.learn,n.componentes.principais)
  semimetrica2 <- SemimetricPCA(curvas$passado.learn,curvas$passado.test,n.componentes.principais)
  
  #        curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,PRED=curvas$passado.test,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values
  #        print(mean((valores.reais - curvas$estimacao)^2))
  
  for (n.vizinhos in vetor.n.vizinhos)
  {
    curvas$estimacao  <- FunopareKnn(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,
                                     PRED=curvas$passado.test,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                     SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values
    
    EscreveResumo(maturidade,curvas,s,valores.reais,truncar,paste("nviz ",n.vizinhos,"-deriv regr ",grau.derivada,sep=""))
    
  }
} # fim da deriv

  

EstimacaoCompletaCorte2 <- function(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s,todas.maturidades,vetor.n.componentes.principais,vetor.grau.derivada){
  
  curvas  <- NULL
  
  for (truncar in vetor.truncar)
    for (maturidade in vetor.maturidade)
      for (s in vetor.s){
        
        valores.reais  <- PreparaValorFuturo(base,percentual.testar,truncar,s,
                                             maturidade)      
        curvas <- PreparaCurvasCorte(base,percentual.testar,truncar,s,
                                     maturidade,todas.maturidades)
                    
        # Estima utilizando a semimétrica PCA
        for (n.componentes.principais in vetor.n.componentes.principais) {
          
          semimetrica1 <- SemimetricPCA(curvas$passado.learn.media,curvas$passado.learn.media,n.componentes.principais)
          semimetrica2 <- SemimetricPCA(curvas$passado.learn.media,curvas$passado.test.media,n.componentes.principais)
          #              curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,PRED=curvas$passado.test.media,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values + curvas$media.test
          #mean((valores.reais - curvas$estimacao)^2)
          curvas$estimacao  <-  NULL
          for (n.vizinhos in vetor.n.vizinhos)
          {
            curvas$estimacao  <- FunopareKnn(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,
                                             PRED=curvas$passado.test.media,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                             SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
              + curvas$media.test
            EscreveResumo(maturidade,curvas,s,valores.reais,truncar,paste("nviz ",n.vizinhos,"-deriv regr ",grau.derivada,sep=""))
          }
          
          
        } # fim da PCA
        
        
        # Estima utilizando a semimétrica derivada
        for (grau.derivada in vetor.grau.derivada){
          
          semimetrica1 <- SemimetricDeriv(curvas$passado.learn.media,curvas$passado.learn.media,grau.derivada,nknot=5,range.grid=c(0,1))
          semimetrica2 <- SemimetricDeriv(curvas$passado.learn.media,curvas$passado.test.media,grau.derivada,nknot=5,range.grid=c(0,1))
          
          curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,PRED=curvas$passado.test.media,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values + curvas$media.test
          for (n.vizinhos in vetor.n.vizinhos)
          {
            curvas$estimacao <- FunopareKnn(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,
                                            PRED=curvas$passado.test.media,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                            SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
              + curvas$media.test
            EscreveResumo(maturidade,curvas,s,valores.reais,truncar,n.vizinhos,paste(passo,",",obs.por.curva,",deriv regr ",grau.derivada,sep=""))
          }
          
        } # fim da deriv
        
      }
}






EstimacaoCompletaPasso <- function(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s,vetor.passo,vetor.obs.por.curva,vetor.n.componentes.principais,vetor.grau.derivada){
  
  curvas  <- NULL
  
  for (truncar in vetor.truncar)
    for (maturidade in vetor.maturidade)
      for (s in vetor.s){
        
        valores.reais  <- PreparaValorFuturo(base,percentual.testar,truncar,s,maturidade)      
        
        for (passo in vetor.passo)
          for (obs.por.curva in vetor.obs.por.curva){
            
            curvas <- PreparaCurvasPasso(base,percentual.testar,truncar,s,maturidade,passo,obs.por.curva)
            
            # Estima utilizando a semimétrica PCA
            for (n.componentes.principais in vetor.n.componentes.principais) {
              
              semimetrica1 <- SemimetricPCA(curvas$passado.learn.media,curvas$passado.learn.media,n.componentes.principais)
              semimetrica2 <- SemimetricPCA(curvas$passado.learn.media,curvas$passado.test.media,n.componentes.principais)
              #              curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,PRED=curvas$passado.test.media,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values + curvas$media.test
              #mean((valores.reais - curvas$estimacao)^2)
              curvas$estimacao  <-  NULL
              for (n.vizinhos in vetor.n.vizinhos)
              {
                curvas$estimacao  <- FunopareKnn(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,
                                                 PRED=curvas$passado.test.media,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                                 SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
                  + curvas$media.test
                EscreveResumo(maturidade,curvas,s,valores.reais,truncar,paste("nviz ",n.vizinhos,"-deriv regr ",grau.derivada,sep=""))
              }
              
              
            } # fim da PCA
            
            
            # Estima utilizando a semimétrica derivada
            for (grau.derivada in vetor.grau.derivada){
              
              semimetrica1 <- SemimetricDeriv(curvas$passado.learn.media,curvas$passado.learn.media,grau.derivada,nknot=5,range.grid=c(0,1))
              semimetrica2 <- SemimetricDeriv(curvas$passado.learn.media,curvas$passado.test.media,grau.derivada,nknot=5,range.grid=c(0,1))
              
              curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,PRED=curvas$passado.test.media,kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values + curvas$media.test
              for (n.vizinhos in vetor.n.vizinhos)
              {
                curvas$estimacao <- FunopareKnn(Response=curvas$futuro.learn.media,CURVES=curvas$passado.learn.media,
                                                PRED=curvas$passado.test.media,neighbour=n.vizinhos,kind.of.kernel="quadratic",
                                                SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
                  + curvas$media.test
                EscreveResumo(maturidade,curvas,s,valores.reais,truncar,n.vizinhos,paste(passo,",",obs.por.curva,",deriv regr ",grau.derivada,sep=""))
              }
              
            } # fim da deriv
            
          }
      }
}


EstimacaoCompletaPasso2 <- function(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s,vetor.passo,vetor.obs.por.curva,vetor.n.componentes.principais,vetor.grau.derivada){
  
  curvas  <- NULL
  
  for (truncar in vetor.truncar)
    for (maturidade in vetor.maturidade)
      for (s in vetor.s){
        
        valores.reais  <- PreparaValorFuturo(base, percentual.testar, truncar,
                                             s, maturidade)      
        for (passo in vetor.passo)
          for (obs.por.curva in vetor.obs.por.curva){
            
            curvas <- PreparaCurvasPasso(base, percentual.testar, truncar, s,
                                         maturidade, passo, obs.por.curva)
            # Estima as distâncias entre as curvas utilizando a semimétrica 
            # baseada em componentes principais
            for (n.componentes.principais in vetor.n.componentes.principais) {
              
              semimetrica1 <- SemimetricPCA(curvas$passado.learn.subtraiult,
                                            curvas$passado.learn.subtraiult,
                                            n.componentes.principais)
              semimetrica2 <- SemimetricPCA(curvas$passado.learn.subtraiult,
                                            curvas$passado.test.subtraiult,
                                            n.componentes.principais)
              
              curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.subtraiult,
                                                 CURVES=curvas$passado.learn.subtraiult,
                                                 PRED=curvas$passado.test.subtraiult,
                                                 kind.of.kernel="quadratic",
                                                 SEMIMETRIC1=semimetrica1,
                                                 SEMIMETRIC2=semimetrica2)$Predicted.values + 
                + curvas$ultimas.obs.test
              
              EscreveResumo(maturidade,curvas,s,valores.reais,nome_arquivos,truncar,paste(passo,",",obs.por.curva,",pca regr ",n.componentes.principais,sep=""))
              
              
            } # fim da PCA
            
            
            # Estima utilizando a semimétrica derivada
            for (grau.derivada in vetor.grau.derivada){
              
              semimetrica1 <- SemimetricDeriv(curvas$passado.learn.subtraiult,
                                              curvas$passado.learn.subtraiult,
                                              grau.derivada,nknot=5,range.grid=c(0,1))
              semimetrica2 <- SemimetricDeriv(curvas$passado.learn.subtraiult,
                                              curvas$passado.test.subtraiult,grau.derivada,nknot=5,range.grid=c(0,1))
              
              curvas$estimacao <- FunopareKnnLcv(Response=curvas$futuro.learn.subtraiult,
                                                 CURVES=curvas$passado.learn.subtraiult,
                                                 PRED=curvas$passado.test.subtraiult,
                                                 kind.of.kernel="quadratic",SEMIMETRIC1=semimetrica1,
                                                 SEMIMETRIC2=semimetrica2)$Predicted.values + 
                + curvas$ultimas.obs.test
              
              #Se nomes de arquivos forem fornecidos, escreve neles; senão, apenas exibe o MSE
              if (!missing(nome_arquivos)) {
                EscreveResumo(maturidade,curvas,s,valores.reais,nome_arquivos,truncar,paste(passo,",",obs.por.curva,",deriv regr ",grau.derivada,sep=""))
              } else {
                mean((valores.reais - curvas$estimacao)^2)
              }
              
            } # fim da deriv
            
          }
      }
}


EstimacaoCompletaRW <- function(base,percentual.testar,vetor.truncar,vetor.maturidade,vetor.s){
  
  curvas  <- NULL
  
  for (truncar in vetor.truncar)
    for (maturidade in vetor.maturidade)
      for (s in vetor.s){
        
        valores.reais  <- PreparaValorFuturo(base,percentual.testar,truncar,s,maturidade)
        curvas <- PreparaCurvasRW(base,percentual.testar,truncar,s,maturidade)
        EscreveResumo(maturidade,curvas,s,valores.reais,truncar)
        
      }
}




ExtraiMediaEUltimaObs <- function (curvas) {
  # Essa função serve para subtrair a média de cada curva e colocar os 
  # valores encontrados na lista curvas
  
  # ARGS
  #    Curva é uma variável criada por alguma função PreparaCurva
  
  curvas$media.learn <- rowMeans(curvas$passado.learn, na.rm = TRUE)
  curvas$media.test <- rowMeans(curvas$passado.test, na.rm = TRUE)
  
  curvas$passado.learn.media <- curvas$passado.learn - curvas$media.learn
  curvas$futuro.learn.media <- curvas$futuro.learn - curvas$media.learn
  curvas$passado.test.media <- curvas$passado.test - curvas$media.test
  
  # Aqui, os valores da última observação são retirados e armazenados em 
  # subtraiult
  obs.por.curva  <- length(curvas$passado.learn[1,])
  curvas$ultimas.obs.learn  <-  curvas$passado.learn[,obs.por.curva]
  curvas$ultimas.obs.test  <-  curvas$passado.test[,obs.por.curva]
  
  curvas$passado.learn.subtraiult  <-  curvas$passado.learn - curvas$ultimas.obs.learn
  curvas$futuro.learn.subtraiult <- curvas$futuro.learn - curvas$ultimas.obs.learn
  curvas$passado.test.subtraiult <- curvas$passado.test - curvas$ultimas.obs.test
  
  return(curvas)
}


FunopareSimples <- function (SEMIMETRIC2, Response, PRED, neighbour) {
  # Faz a estimação Kernel  
  ordem <- apply(X=SEMIMETRIC2,MARGIN=2,FUN=order)
  
  h <- semimetrica2[apply(X=semimetrica2,MARGIN=2,FUN=order)[knn, ]]
  sum(semimetrica2[,30] <= h[30])
  
  estimacao  <- rep(0,ncol(semimetrica2))
  for (i in ncol(semimetrica2)){
    indices.menor.distancia <- which(order(semimetrica2[, i]) < knn)
    estimacao[i] <- Response[indices.menor.distancia]*semimetrica2[indices.menor.distancia, i] / sum(semimetrica2[indices.menor.distancia, i] / h.knn)
  }
  return(estimacao)
}


FunopareKnn <- function(Response, CURVES, PRED, neighbour, kind.of.kernel = "quadratic", SEMIMETRIC1, SEMIMETRIC2)
{
  ################################################################
  # Função modificada para agilidade quando calculada para múltiplos parâmetros.
  # As duas funções semimétrica já devem ser fornecidas, sendo SEMIMETRIC1 as distancias
  # d(learning,learning) e SEMIMETRIC2 as distancias de d(learning,testing)
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A bandwidth corresponding to number of neighbours has to be given.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "neighbour" number of neighbours fixed for computing the
  #                functional kernel estimator.
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "kNN" value of the current argument "neighbour".
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  kernel <- get(kind.of.kernel)
  p1 <- ncol(SEMIMETRIC1)
  n1 <- nrow(SEMIMETRIC1)
  if(neighbour >= n1)
    neighbour <- n1
  bandwidth.knn1 <- 0
  for(j in 1:p1) {
    Sem <- SEMIMETRIC1[, j]
    knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
    bandwidth.knn1[j] <- 0.5 * sum(knn.to.band)
  }
  KERNEL1 <- kernel(t(t(SEMIMETRIC1)/bandwidth.knn1))
  KERNEL1[KERNEL1 < 0] <- 0
  KERNEL1[KERNEL1 > 1] <- 0
  diag(KERNEL1) <- 0
  RESPKERNEL1 <- KERNEL1 * Response
  Denom1 <- apply(KERNEL1, 2, sum)
  Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    p2 <- ncol(SEMIMETRIC2)
    bandwidth.knn2 <- 0
    for(j in 1:p2) {
      Sem <- SEMIMETRIC2[, j]
      knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)
                                    ]]
      bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
    }
    KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
    KERNEL2[KERNEL2 < 0] <- 0
    KERNEL2[KERNEL2 > 1] <- 0
    Denom2 <- apply(KERNEL2, 2, sum)
    RESPKERNEL2 <- KERNEL2 * Response
    Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, knn = neighbour, 
                Mse = Mse.estimated))
  }
  else {
    return(list(Estimated.values = Response.estimated, knn = 
                  neighbour, Mse = Mse.estimated))
  }
}

FunopareKnnGcv <- function(Response, CURVES, PRED, kind.of.kernel = "quadratic", SEMIMETRIC1, SEMIMETRIC2)
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A global bandwidth (i.e. a number of neighbours) is selected by a 
  # cross-validation procedure.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Bandwidths" vector containing the global data-driven bandwidths  
  #                 for each curve in the matrix "CURVES"
  #    "knearest.opt" optimal number of neighbours
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  step <- ceiling(n1/100)
  if(step == 0)
    step <- 1
  Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  kmax <- max(Knearest)  
  # the vector Knearest contains the sequence of the 
  # k-nearest neighbours used for computing the optimal bandwidth
  Response.estimated <- 0
  Bandwidth.opt <- 0
  HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
  BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]  
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)	
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]	
    # BANDWIDTH[i, l-1] contains (dq(X_{j_l},X_i) + 
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/BANDWIDTH[i,  ]
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0
    Ind.curves <- Norm.order[2:(kmax + 1)]
    Ind.resp <- Response[Ind.curves]
    YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
    HAT.RESP[i,  ] <- apply(YMAT[Knearest,  ] * KMAT[Knearest,  ], 
                            1, sum)/apply(KMAT[Knearest,  ], 1, sum)
  }
  CRITERIUM <- (HAT.RESP - Response)^2
  Criterium <- apply(CRITERIUM, 2, sum)
  index.opt <- order(Criterium)[1]
  Response.estimated <- HAT.RESP[, index.opt]
  knearest.opt <- Knearest[index.opt]
  Bandwidth.opt <- BANDWIDTH[, knearest.opt]
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(KERNEL, 2, sum)
    RESPKERNEL <- KERNEL * Response
    Response.predicted <- apply(RESPKERNEL, 2, sum)/Denom
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, Bandwidths = 
                  Bandwidth.opt, knearest.opt = knearest.opt, Mse = 
                  Mse.estimated))
  }else {
    return(list(Estimated.values = Response.estimated, Bandwidths
                = Bandwidth.opt, knearest.opt = knearest.opt, Mse = 
                  Mse.estimated))
  }
}

FunopareKnnLcv <- function(Response, CURVES, PRED, kind.of.kernel = "quadratic", SEMIMETRIC1, SEMIMETRIC2){
  ################################################################
  # Função modificada para agilidade quando calculada para múltiplos parâmetros.
  # As duas funções semimétrica já devem ser fornecidas, sendo SEMIMETRIC1 as distancias
  # d(learning,learning) e SEMIMETRIC2 as distancias de d(learning,testing)
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A local bandwidth (i.e. local number of neighbours) is selected 
  # by a cross-validation procedure.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Bandwidths" vector containing the local data-driven bandwidths
  #                 for each curve in the matrix "CURVES"
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  step <- ceiling(n1/100)
  if(step == 0)
    step <- 1
  Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  kmax <- max(Knearest)
  # the vector Knearest contains the sequence of the 
  # k-nearest neighbours used for computing the optimal bandwidth
  Response.estimated <- 0
  Bandwidth.opt <- 0
  Knn1 <- 0
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]
    # Bandwidth[l-1] contains (dq(X_{j_l},X_i) + 
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    Bandwidth <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/Bandwidth
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0
    Ind.curves <- Norm.order[2:(kmax + 1)]
    Ind.resp <- Response[Ind.curves]
    YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
    Hat.resp <- apply(YMAT[Knearest,  ] * KMAT[Knearest,  ], 1, sum
    )/apply(KMAT[Knearest,  ], 1, sum)
    Criterium <- abs(Hat.resp - Response[i])
    index <- order(Criterium)[1]
    Knn1[i] <- Knearest[index]
    Response.estimated[i] <- Hat.resp[index]
    Bandwidth.opt[i] <- Bandwidth[index]
  }
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Sm2k.ord <- order(SEMIMETRIC2[, k])
      knn <- Knn1[Sm2k.ord[1]]
      Bandwidth2[k] <- sum(sort(Sm2k)[knn:(knn+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(as.matrix(KERNEL), 2, sum)
    RESPKERNEL <- KERNEL * Response
    Response.predicted <- apply(as.matrix(RESPKERNEL), 2, sum)/
      Denom
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, Bandwidths = 
                  Bandwidth.opt, Mse = Mse.estimated))
  }
  else {
    return(list(Estimated.values = Response.estimated, Bandwidths
                = Bandwidth.opt, Mse = Mse.estimated))
  }
}


# GraficoImagem <- function(df, ...){
# # A partir de um data.frame cria um gráfico de calor
# # O data.frame deve possuir como primeiros argumentos seus parâmetros, e
# #   na sequência apenas os intervalos de análise
#   nome.cols <- names(df)
#   colunas.intervalos <- grep(pattern="[0-9]*-[0-9]*",nome.cols)  
#   which.min(apply(df[, colunas.intervalos],1,sum))
#   
#   for(i in (setdiff(1:ncol(df),colunas.intervalos))
#     if (is.factor(df[, i]))
#       df[, i] <- as.integer(df[, i])
#   
#   str(df)
# }

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

# Multiplot.ts(matrix(rep(1:20,each=20),nrow=20))

Multiplot.ts <- function(data, ... , range){
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
  plot(data[, range[1]],col=1,ylim=c(min(data[, range]), 
                                     max(data[, range])), type="l",  ...)
  for (i in range[-1]) {
    lines(data[, i],col=i, lty=((i - 1) %/% 8 + 1), ...)   
  }
    
  
}


Multiplotdiv.ts <- function(data, div){
  
  #### TO FINISH ####
  
  # This function plots many time series in the same graph
  #
  # ARGS  
  #    Data is a matrix with each serie in its lines
  #    Div is an optional vector with the indexes of the time series
  #       to be exibited. If it is absent, then every line is ploted
  
  nlinhas  <- nrow(data)
  j <- 1
  while (j <= nlinhas) {
  
    plot.ts(data[range[1],],col=1,ylim=c(min(data[range,]),max(data[range,])))
    for (i in j:(j+div))
      lines.ts(data[i,],col=i)  
  }
}

py.range <- function(range){
  # implements python function range, but starting on 1 instead of 0
  return(range[1]:range[2])
}

PreparaCurvasCorte  <- function(base,percentual.testar,intervalo,s,maturidade,todas.maturidades = NULL,retirar = NULL) {
########################################################################
# Essa função retorna uma variável com as curvas necessárias para fazer a estimação.
# 
# ARGS
#   base: uma base de dados; cada linha representa uma série temporal
#   percentual.testar: o percentual da amostra que deve ser utilizado para testar 
#         a capacidade de previsão do modelo
#   intervalo: o intervalo de análise da série temporal
#   s: o horizonte de previsão
#   maturidade: índice da maturidade que se deseja prever
#   todas.maturidades: indica as maturidades que se quer utilizar para compor as curvas
#   retirar: vetor com as quantidades que devem ser retiradas do valor base para
#               cálculo da semimétrica (o índice dos intervalos de retirar deve
#               ser os mesmos da base)
#####################################################################################
  if(missing(retirar))
    retirar <- rep(0,dim(base)[1])
  if(missing(todas.maturidades))
    todas.maturidades <- 1:dim(base)[2]
  if ((dim(base)[2]) < max(todas.maturidades) ||  min(todas.maturidades) < 1)  
    stop("Você forneceu um conjunto inválido de maturidades na variável todas.maturidades")  
  if (intervalo[1] < 1 || intervalo[2] > dim(base)[1])
    stop("o intervalo para truncar a série temporal é inválido")
  if (0 >= percentual.testar || percentual.testar >= 1)
    stop("O percentual de teste fornecido é inválido")
  tamanho.serie.temporal  = intervalo[2] - intervalo[1] + 1  
  #learning[1] contém o indice da primeira curva para aprendizado e learning[2] o índice da última
  learning <- c(1,trunc((1-percentual.testar)*tamanho.serie.temporal)-s) + intervalo[1] - 1
  #testing é igual ao learning, mas para as curvas usadas para previsão
  testing <- c(trunc((1-percentual.testar)*tamanho.serie.temporal),tamanho.serie.temporal-s) + intervalo[1] - 1
  if(learning[1] >= learning[2] || testing[1] >= testing[2])
    stop("Você forneceu um valor inválido de horizonte de previsão")  
  curvas <- NULL
  curvas$retirar.learn <- retirar[learning[1]:learning[2]]
  curvas$retirar.test <- retirar[testing[1]:testing[2]]
  curvas$passado.learn <- as.matrix(base[learning[1]:learning[2], ] - retirar[learning[1]:learning[2]])
  curvas$passado.test <- as.matrix(base[testing[1]:testing[2], ] - retirar[testing[1]:testing[2]])
  curvas$futuro.learn <- as.matrix(base[learning[1]:learning[2]+s,maturidade] - retirar[learning[1]:learning[2]])
  curvas$estimacao  <- NULL
  curvas$tipo  <- "corte"
  class(curvas) <- "fdaCorte"
  return(curvas)
}


predict.fdaCorte <- function(object, semimetricas, n.vizinhos = NULL, cv = "GCV") {
###########################################################################
# Função que faz previsões para objetos da classe fdaCorte
# ARGS
#   object: um objeto da classe fda.Corte
#   n.vizinhos: valor opcional para determinar ou não a quantidade de vizinhos
#                            mais próximos para se fazer a estimação
#   cv: como é feito cross-validation para determinar o número de vizinhos ótimo knn.
#             aceita dois tipos:  "GCV" (default) - cross-validation global
#                                 "LCV"- cross-validation local
# RETORNA
#   vetor com os valores futuros previstos
#############################################################################
  if(class(semimetricas) != "semimetricas")
    stop("Você não forneceu curvas de semimétrica válidas")
  semimetrica1 = semimetricas[[1]]; semimetrica2 = semimetricas[[2]]
  if (missing(n.vizinhos)){
    if (cv == "GCV") {
      estimacao  <- FunopareKnnGcv(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,
                                PRED=curvas$passado.test,kind.of.kernel="quadratic",
                                SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
                    curvas$retirar.test
    } else {
      if (cv != "LCV")
        stop("O tipo de cross-validation fornecido é inválido")
      estimacao  <- FunopareKnnLcv(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,
                                   PRED=curvas$passado.test,kind.of.kernel="quadratic",
                                   SEMIMETRIC1=semimetrica1,SEMIMETRIC2=semimetrica2)$Predicted.values +
        curvas$retirar.test
    }
  } else {
    estimacao  <- FunopareKnn(Response=curvas$futuro.learn,CURVES=curvas$passado.learn,
                              PRED=curvas$passado.test,neighbour=n.vizinhos,
                              kind.of.kernel="quadratic", SEMIMETRIC1=semimetrica1,
                              SEMIMETRIC2=semimetrica2)$Predicted.values + 
                  curvas$retirar.test
  }
  return(estimacao)
}


PreparaCurvasRW  <- function(base,percentual.testar,truncar,s,maturidade) {
  # Essa função retorna uma variável com os valores previstos pelo método RandomWalk
  
  curvas <- NULL
  
  if (truncar[1] < 1 || truncar[2] > dim(base)[1])
    stop("o intervalo para truncar a série temporal é inválido")
  
  if (0 >= percentual.testar || percentual.testar >= 1)
    stop("O percentual de teste fornecido é inválido")
    
  tamanho.serie.temporal  = truncar[2]-truncar[1]+1
    
  #testing é igual ao learning, mas para as curvas usadas para previsão
  testing <- c(trunc((1-percentual.testar)*tamanho.serie.temporal) + truncar[1]-1,tamanho.serie.temporal - s + truncar[1] - 1)
  
  if(testing[1] >= testing[2])
    stop("Você forneceu um valor inválido de horizonte de previsão")
  
  curvas$estimacao <- base[testing[1]:testing[2],maturidade+1] 
  curvas$tipo  <- "rw"
  
  return(curvas)
  
}


PreparaCurvasPasso  <- function(base,percentual.testar,truncar,s,maturidade,
                                passo,obs.por.curva) {
  
  curvas  <- NULL
  
  
  if(passo < 1)
    stop(paste("Você forneceu um passo inválido"))
    
  if (truncar[1] < 1 || truncar[2] > dim(base)[1])
    stop("o intervalo para truncar a série temporal é inválido")
  
  if (0 >= percentual.testar || percentual.testar >= 1)
    stop("O percentual de teste fornecido é inválido")
  
  z <- ts(data=base[truncar[1]:truncar[2],maturidade])
  tamanho.serie.temporal  = length(z)
  
  learning <- c(1, trunc((1-percentual.testar)*tamanho.serie.temporal) - 
                  obs.por.curva-s+1)
  
  testing <- c(trunc((1-percentual.testar)*tamanho.serie.temporal) + 
                 1 - obs.por.curva,tamanho.serie.temporal-s-obs.por.curva+1)
  
  curvas$passado.learn <- NULL
  curvas$passado.test <- NULL
  curvas$futuro.learn <- NULL
  
  i = learning[2] # o índice i representa a primeira observação das curvas que serão armazenadas
  # laço vai preenchendo as curvas de tras para frente
  while (i >= learning[1]) 
  {
    curvas$passado.learn <- cbind(z[i:(i+obs.por.curva-1)],curvas$passado.learn)  
    curvas$futuro.learn <- c(z[i+obs.por.curva-1+s],curvas$futuro.learn)
    i = i - passo 
  }
  curvas$passado.learn <- t(curvas$passado.learn)
  
  # o índice j é índice da matriz y, enquanto i é índice da série temporal z
  for (i in testing[1]:testing[2])
  {
    curvas$passado.test <- cbind(curvas$passado.test,z[i:(i+obs.por.curva-1)]) 
  }
  curvas$passado.test <- t(curvas$passado.test)
  curvas$tipo  <- "passo"
  
  curvas <- ExtraiMediaEUltimaObs(curvas)
    
  return(curvas)
  
}


PreparaValorFuturo  <-function (base,percentual.testar,truncar,s,maturidade) {
# Esta função retorna os valores verdadeiros, para comparar aos valores
# estimados pelos procedimentos
  if (truncar[1] < 1 || truncar[2] > dim(base)[1])
    stop("o intervalo para truncar a série temporal é inválido")
  if (0 >= percentual.testar || percentual.testar >= 1)
    stop("O percentual de teste fornecido é inválido")
  tamanho.serie.temporal  = truncar[2]-truncar[1]+1
  #testing é igual ao learning, mas para as curvas usadas para previsão
  testing <- c(trunc((1-percentual.testar)*tamanho.serie.temporal) + 
                 truncar[1]-1,tamanho.serie.temporal - s + truncar[1] - 1)
  if(testing[1] >= testing[2])
    stop("Você forneceu um valor inválido de horizonte de previsão")  
  valores.reais <- base[testing[1]:testing[2]+s,maturidade] 
  return(valores.reais)
}


SemimetricasClasse <- function(curvas, ..., tipo = "pca") {
################################################################  
# Prepara um objeto da classe "semimetricas", para ser usado para previsão
# ARGS
#    curvas: um objeto do tipo curvas
#    tipo: o tipo das semimetricas - "pca" (default) ou "deriv"
#     ... : argumentos restantes para o cálculo das semimétricas
#       - pca: q: número de componentes principais. 
#                 valores default: q = 5
#       - deriv:  q: ordem da derivação
#                 nknot: número de nós para cálculo do spline 
#                 range.grid: vetor de tamanho 2 contendo o intervalo
#                             no qual as curvas serão discretizadas
#                 valores default: q = 0, nknot = 5, range.grid = c(0,10)
# RETORNA
#   uma lista com duas matrizes:
#   semimetrica1: matriz contendo as distâncias entre cada duas curvas de treino
#   semimetrica2: matriz contendo as distâncias entre cada curva de teste com
#                           cada curva de treino
################################################################  
  if (tipo == "pca"){
    semimetrica1 <- SemimetricPCA(curvas$passado.learn,
                                  curvas$passado.learn, ...)
    semimetrica2 <- SemimetricPCA(curvas$passado.learn,
                                  curvas$passado.test, ...)  
  } else {
    if (tipo != "deriv") 
      stop("Você forneceu um tipo inválido de semimétrica")
    semimetrica1 <- SemimetricDeriv(curvas$passado.learn,
                                    curvas$passado.learn, ...)
    semimetrica2 <- SemimetricDeriv(curvas$passado.learn,
                                    curvas$passado.test, ...)  
  }
  semimetricas <- list(semimetrica1,semimetrica2)
  class(semimetricas) <- "semimetricas"
  return(semimetricas)
}


SemimetricPCA <- function(DATA1, DATA2, q = 5){
  ###############################################################
  # Computes between curves a pca-type semimetric based on the
  # functional principal components analysis method.
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "q" the retained number of principal components
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  qmax <- ncol(DATA1)
  if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
  n <- nrow(DATA1)
  COVARIANCE <- t(DATA1) %*% DATA1/n
  EIGENVECTORS <- eigen(COVARIANCE, sym = T)$vectors[, 1:q]
  COMPONENT1 <- DATA1 %*% EIGENVECTORS
  if(twodatasets) {
    COMPONENT2 <- DATA2 %*% EIGENVECTORS
  }
  else {
    COMPONENT2 <- COMPONENT1
  }
  SEMIMETRIC <- 0
  for(qq in 1:q)
    SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, qq], COMPONENT2[, 
                                                                  qq], "-")^2
  return(sqrt(SEMIMETRIC))
}


SemimetricDeriv <- function(DATA1, DATA2, q = 0, nknot = 5, range.grid = c(0,10)){
  ###############################################################
  # Computes a semimetric between curves based on their derivatives.
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "q" order of derivation
  #    "nknot" number of interior knots (needed for defining the B-spline basis)
  #    "range.grid" vector of length 2 containing the range of the grid at 
  #                 which the curve are evaluated (i.e. range of the 
  #                 discretization)
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  library(splines)
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  #####################################################################
  # B-spline approximation of the curves containing in DATASET :
  # -----------------------------------------------------------
  # "knot" and "x" allow to define the B-spline basis
  # "coef.mat1[, i]" corresponds to the B-spline expansion
  # of the discretized curve contained in DATASET[i, ]. 
  # The B-spline approximation of the curve contained in "DATA1[i, ]" 
  # is given by "Bspline %*% coef.mat1[, i]"
  #####################################################################
  p <- ncol(DATA1)
  a <- range.grid[1]
  b <- range.grid[2]
  x <- seq(a, b, length = p)
  order.Bspline <- q + 3
  nknotmax <- (p - order.Bspline - 1)%/%2
  if(nknot > nknotmax){
    stop(paste("give a number nknot smaller than ",nknotmax, " for avoiding ill-conditioned matrix"))
  }
  Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
  delta <- sort(c(rep(c(a, b), order.Bspline), Knot))
  Bspline <- splineDesign(delta, x, order.Bspline)
  Cmat <- crossprod(Bspline)
  Dmat1 <- crossprod(Bspline, t(DATA1))
  coef.mat1 <- symsolve(Cmat, Dmat1)
  #######################################################################
  # Numerical integration by the Gauss method :
  # -------------1------------------------------
  # The objects ending by "gauss" allow us to compute numerically  
  # integrals by means the "Gauss method" (lx.gauss=6 ==> the computation 
  # of the integral is exact for polynom of degree less or equal to 11).
  #######################################################################
  point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 
                   0.2386191861, 0.6612093865, 0.9324695142)
  weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
  x.gauss <- 0.5 * ((b + a) + (b - a) * point.gauss)
  lx.gauss <- length(x.gauss)
  Bspline.deriv <- splineDesign(delta, x.gauss, order.Bspline, rep(q, lx.gauss))
  H <- t(Bspline.deriv) %*% (Bspline.deriv * (weight.gauss * 0.5 * (b - a)))
  eigH <- eigen(H, sym = T)
  eigH$values[eigH$values < 0] <- 0
  Hhalf <- t(eigH$vectors %*% (t(eigH$vectors) * sqrt(eigH$values)))
  COEF1 <- t(Hhalf %*% coef.mat1)
  if(twodatasets){
    Dmat2 <- crossprod(Bspline, t(DATA2))
    coef.mat2 <- symsolve(Cmat, Dmat2)
    COEF2 <- t(Hhalf %*% coef.mat2)
  } else {
    COEF2 <- COEF1
  }
  SEMIMETRIC <- 0
  nbasis <- nrow(H)
  for(f in 1:nbasis)
    SEMIMETRIC <- SEMIMETRIC + outer(COEF1[, f], COEF2[, f], "-")^2
  return(sqrt(SEMIMETRIC))
}
