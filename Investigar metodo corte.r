# Testes ------------------------------------------------------------------
i <- 1

curvas <- PreparaCurvasCorte(taxas.juro,percentual.testar,c(1,2500),s=61,maturidade=2,todas.maturidades,betas[, 1])
valor.real <- PreparaValorFuturo(taxas.juro,percentual.testar,c(1,2500),s=61,maturidade=2)
# curvas1 <- PreparaCurvasCorte(taxas.juro,percentual.testar,c(1001,2000),s=121,maturidade=2,todas.maturidades,taxas.juro[, 1])
# str(curvas1)
# View(curvas$passado.learn)
# View(curvas$passado.test)
Distancias <- SemimetricPCA(curvas$passado.learn,curvas$passado.test,7)
#View(Distancias)
# Multiplot.ts(t(curvas$passado.learn))
# Multiplot.ts(t(curvas1$passado.learn))
dist <- 3

# tmp <- as.data.frame(cbind(Distancias=Distancias[which(Distancias[, i] < dist), i],
#              Resposta = curvas$futuro.learn[which(Distancias[, i] < dist)],
#              Diff.resp = abs(curvas$futuro.learn[which(Distancias[, i] < dist)] - 
#                                curvas$futuro.learn[i]),
#              Indice = which(Distancias[, i] < dist)))


# compara a relação entre as respostas dos testes com a das respostas
tmp2 <- as.data.frame(cbind(Distancias=Distancias[which(Distancias[, i] < dist), i],
                            Resposta = curvas$futuro.learn[which(Distancias[, i] < dist)],
                            Diff.resp = abs(curvas$futuro.learn[which(Distancias[, i] < dist)] - 
                                              valor.real[i] + curvas$retirar.learn[which(Distancias[, i] < dist)]) ,
                            Indice = which(Distancias[, i] < dist)))

View(tmp2)

qplot(Distancias,Diff.resp, data= tmp2)

summary(lm(Diff.resp ~ Distancias,tmp2))


# 
# View(tmp[order(tmp$Distancias), ])
# Multiplot.ts(t(curvas$passado.learn[which(Distancias[, i] < dist), ]))
# 
