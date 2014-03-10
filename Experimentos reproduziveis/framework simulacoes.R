
# Bibliotecas para carregar -----------------------------------------------

# Carrega as funções para visualização e cálculo dos intervalos O pacote abaixo
# está disponível para download no endereço abaixo:
# https://www.dropbox.com/s/0v71p5soytpfqqm/biblioteca%20simulacoes.R
source("biblioteca simulacoes.R") 


# outras bibliotecas (...)



# Parâmetros para realizar os experimentos --------------------------------

# intervalos: uma lista com os intervalos nos quais se fará a análise
vetor.intervalos <- list(c(1,800),c(201,1000),c(401,1200))

# vetor.maturidade: contém os índices das maturidades que serão estimadas 
# e previstas
vetor.maturidade = c(1:4,6,8,12,17)

# vetor.horizontes: contém todos os horizontes de previsão para testar
vetor.horizontes <- c(1,2,3)


# A variável percentual.testar determina quanto do intervalo total será usado
# para testar; 0.5 = 50% É utilizado para definir as variáveis intervalo.passado
# e intervalo.futuro. Essas duas variáveis são necessárias para dar entrada nas
# funções de previsão. intervalo.passado é um vetor de dimensão 2 que indica o
# início e o fim do intervalo que será utilizado como treino. Já
# intervalo.futuro também é um vetor de dimensão 2, indicando o início e o fim
# do intervalo que a previsão será feita. Exemplo: sendo o valor das duas
# variáveis, respectivamente, c(1,1000) e c(1200,1400), a estimação é realizada
# no intervalo 1:1000 e 201 valores serão previstos, para o intervalo 1200:1400.
# Finalmente, o que a variável percentual testar faz, portanto, é definir o
# tamanho relativo entre o intervalo.passado e intervalo.futuro. Consulte também
# a documentação das funções IntervaloPassado e IntervaloFuturo.
percentual.testar <- 0.2


# Simulações --------------------------------------------------------------


simulacoes <- vector("list", 0) # RESETA os dados de simulacoes

for (horizonte in vetor.horizontes) {
  for (j in 1:length(vetor.intervalos)) {
    
    intervalo.total <- vetor.intervalos[[j]] 
    intervalo.passado <- IntervaloPassado(intervalo.total,percentual.testar) # ler a documentação dessas funções
    intervalo.futuro <- IntervaloFuturo(intervalo.total,percentual.testar, horizonte)

    ## Métodos que prevêem simultaneamente para todas as maturidades, como o
    ## Diebold-li ou VAR(1), ficam aqui. Dentro do próximo laço (da maturidade)
    ## que os valores previstos serão gravados ##

    # Método A #
    previsao.metodoA <- FuncaoMetodoA(... ,intervalo, horizonte)
    
    # Método B # 
    previsao.metodoB <- FuncaoMetodoB(... ,intervalo, horizonte)
    
    # (...) #
    
    for (i in 1:length(vetor.maturidade)) {
      maturidade <- vetor.maturidade[i]

      # Método 1 #
      previsao.metodo1 <- FuncaoMetodo1(... , maturidade, intervalo, horizonte)
      
      # Método 2 # 
      previsao.metodo2 <- FuncaoMetodo2(... , maturidade, intervalo, horizonte)
      
      # (...) #
      
      
      ## Valores reais ##
      valores.reais <- taxas.juro[intervalo.futuro[1]:intervalo.futuro[2], maturidade]
      
      ## Grava no arquivo ##
      
      # Para as previsões do laço de fora, deve-se retirar apenas a maturidade desejada
      simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["metodoA"]]  <- previsao.metodoA[, maturidade]
      simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["metodoB"]] <- previsao.metodoB[, maturidade]
      
      # Para as outras, a estimação toda pode ser gravada
      simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["metodo1"]] <- previsao.metodo1
      simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["metodo2"]] <- previsao.metodo2
      simulacoes[[Nome.vetor(intervalo.total)]][[as.character(horizonte)]][[as.character(maturidade)]][["valores.reais"]] <- valores.reais
    }
  }
}

# salva o objeto simulações no lugar desejado. Aconselha-se salvar com
# frequência para não correr o risco de sobrescrever os experimentos e perder
# todo o trabalho
save(simulacoes,file="simulacoes_reduzida.dad")



# Análise dos resultados --------------------------------------------------


# Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
# de maturidade x horizonte x intervalo
tabela <- Tabelao(simulacoes,vetor.maturidade,c(5,25,50),list(c(1,800),c(201,1000)))
heatmap(tabela)

# Exibe 3 gráficos: I) séries previstas por cada método e a realizada; II) erro
# quadrático de cada série prevista, por unidade de tempo; III) erro quadrático
# médio de cada série, até a unidade de tempo t. 
PlotarSeries(simulacoes,maturidade=2,horizonte=5,intervalo=c(1,800))