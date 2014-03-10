
# Análise dos arquivos ----------------------------------------------------

mat <- "mat3"
hor <- "hor32"
i <- 2


library(plyr)
library(reshape2)

diretorio.analise <- "./Resultados/14-out/"

# setwd("/media/marcelo/ee2881fc-31fa-490d-82f7-fc4ed9292280/home/ruas/Dropbox/R/NPFDA/Resultados/6-nov")

# Linhas contém o nome dos arquivos no diretódio
linhas <- as.vector(list.files(diretorio.analise))

# O trecho abaixo procura no diretório de análise quais são 
ocorrencias <- regexpr("^mat[0-9]*",linhas)
(maturidades.existentes <- levels(as.factor(regmatches(linhas,ocorrencias))))
ocorrencias <- regexpr("hor[0-9]*",linhas)
(horizontes.existentes <- levels(as.factor(regmatches(linhas,ocorrencias))))
ocorrencias <- regexpr("[0-9]*-[0-9]*",linhas)
(intervalos.existentes <- levels(as.factor(regmatches(linhas,ocorrencias))))



# O laço abaixo seleciona as combinações de maturidades e horizontes
for (mat in maturidades.existentes){
  for (hor in horizontes.existentes) {
    arquivos.matur.horiz <- grep(paste(mat,"_",hor,sep=""),list.files(diretorio.analise))
    
    # cria linhas com valores nulos para serem preenchidos após -- INUTILIZADO
#     linha.passo <- rep(NA,8+length(intervalos.existentes))
#     names(linha.passo) <- c("Maturidade","Horizonte","KNN","PASSO", "OBS/CURVA", "TIPO SEMIMETRICA",intervalos.existentes,"Vitorias sobre RW","Total")
#     linha.passo[1:2] <- c(mat,hor)
#     
#     linha.corte <- rep(NA,6+length(intervalos.existentes))
#     names(linha.passo) <- c("Horizonte","Maturidade","KNN", "TIPO SEMIMETRICA",intervalos.existentes,"Vitorias sobre RW","Total")
#     
#     linha.rw <- rep(NA,3+length(intervalos.existentes))
#     names(linha.passo) <- c("Horizonte","Maturidade",intervalos.existentes,"Total")
  
    # cria arquivo com os resultados para uma dada maturidade e horizonte
    arquivo.para.escrever <- paste0(diretorio.analise,"Relatorios/Relatorio_",mat,"_",hor,".txt")
    
    
    # cria a tabela para os resultados a partir do primeiro arquivo
    nome.arquivo <- list.files(diretorio.analise)[arquivos.matur.horiz[1]]
    arquivo <- read.csv(paste0(diretorio.analise,nome.arquivo), header=FALSE)
    # extrai o intervalo (ex: 1201-2200) do nome do arquivo
    intervalo <- unlist(strsplit(unlist(strsplit(nome.arquivo, "[.]"))[1],"_"))[3]
    
    tabela.passo <- arquivo[arquivo[, 2] == "passo", c(3:6,1)]
    names(tabela.passo) <- c("KNN","PASSO", "OBS/CURVA", "TIPO SEMIMETRICA",intervalo)
    
    tabela.corte <- arquivo[arquivo[, 2] == "corte", c(3:4,1)]
    names(tabela.corte) <- c("KNN", "TIPO SEMIMETRICA",intervalo)
    
    tabela.rw <- arquivo[arquivo[, 2] == "rw", c(1)]
    names(tabela.rw) <- c(intervalo)
    
    for (i in arquivos.matur.horiz[-1]) {
      
            
      #   seleciona o i-ésimo arquivo com a maturidade=mat e horizonte=hor e executa 
      # os msms comandos de antes, dessa vez para o i-ésimo intervalo
      nome.arquivo <- list.files(diretorio.analise)[i]
      arquivo <- read.csv(paste0(diretorio.analise,nome.arquivo), header=FALSE)
      intervalo <- unlist(strsplit(unlist(strsplit(nome.arquivo, "[.]"))[1],"_"))[3]
      
      tabela.passo.i <- arquivo[arquivo[, 2] == "passo", c(3:6,1)]
      names(tabela.passo.i) <- c("KNN","PASSO", "OBS/CURVA", "TIPO SEMIMETRICA",intervalo)
      tabela.corte.i <- arquivo[arquivo[, 2] == "corte", c(3:4,1)]
      names(tabela.corte.i) <- c("KNN", "TIPO SEMIMETRICA",intervalo)
      tabela.rw.i <- arquivo[arquivo[, 2] == "rw", c(1)]
      names(tabela.rw.i) <- c(intervalo)
      
      # une o i-ésimo intervalo ao restante
      tabela.passo <- merge(tabela.passo, tabela.passo.i)
      tabela.corte <- merge(tabela.corte, tabela.corte.i)
      tabela.rw <- c(tabela.rw, tabela.rw.i)
      
            
    }
    tabela.passo
    colunas.eqm.passo <- 1:length(tabela.rw) + 4
    colunas.eqm.corte <- 1:length(tabela.rw) + 2

    # vitórias sobre o RW
    rowSums(tabela.passo[, colunas.eqm.passo] < tabela.rw)
    rowSums(tabela.corte[, colunas.eqm.corte] < tabela.rw)
    
    # Soma dos EQM
    rowSums(tabela.passo[, colunas.eqm.passo]) < sum(tabela.rw)
    rowSums(tabela.corte[, colunas.eqm.corte]) < sum(tabela.rw)
    
    
    ### Analisa os melhores parâmetros para o método passo
    parametros.passo <- arquivo[arquivo[, 2] == "passo", c(3:6,1)]
    
    
    ### Analisa os melhores parâmetros para o método corte
    
    
    ### Imprime os resultados do método rw
          
    # Escreve no arquivo
    cat(linha,file=arquivo.para.escrever,sep="\n",append=TRUE)
    
    
  }
}


write.csv(tabela.passo,"tabela resultados.csv")
