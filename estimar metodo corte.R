source("biblioteca-procedimento-npfda.r")

# intervalo.treino <- 1:2500
# intervalo.teste <- 2501:5000
# intervalo.validacao <- 5001:6530

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
# 3  6  9	12	15	18	21	24	30	36	48	60	72	84	96	108	120	144	180	240
# maturidades <- c(0.25 , 0.5,  0.75,  1,	1.25,	1.5,	1.75,	2,	2.5,	3,	4,	5,	6,	7,	8,	9,	10,	12,	15, 20)
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,	1.5,	1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10,	12, 15, 20)
percentual.testar <- 0.2
intervalo <- 1:6535
taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                        [intervalo ,2:(length(tempo.maturidades)+1)])
datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
                              [intervalo, 1]), format = "%Y-%m-%d")
retirar <- taxas.juro[, 4]

calor <- 400

curvas <- PreparaCurvasCorte(taxas.juro,percentual.testar,intervalo.total,horizonte,maturidade=2,retirar=retirar)
semimetricas <- SemimetricasClasse(curvas,q=3)
taxas.previstas.corte <- predict(curvas,semimetricas,n.vizinhos=30)
