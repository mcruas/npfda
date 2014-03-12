# Testes ------------------------------------------------------------------
maturidade <- 2
intervalo.total <- c(1,2500)

retirar <- taxas.juro[, 10]
curvas <- PreparaCurvasCorte(taxas.juro,intervalo.passado,intervalo.futuro,maturidade=maturidade, retirar = retirar)
semimetricas <- SemimetricasClasse(curvas,q=1,tipo="deriv")
taxas.previstas.corte <- predict(curvas,semimetricas)

intervalo.passado <- IntervaloPassado(intervalo.total,percentual.testar)
intervalo.futuro <- IntervaloFuturo(intervalo.total,percentual.testar, horizonte) 
valores.reais <- taxas.juro[py.range(intervalo.futuro), ]
betas.previstos <- DieboldLi.PreveBetas(betas, intervalo.passado,
                                        intervalo.futuro)
taxas.previstas.diebold <- DieboldLi.BetasParaTaxas(
  betas.previstos, tempo.maturidades, parametro)

taxas.previstas.rw <- taxas.juro[py.range(intervalo.futuro) - horizonte, maturidade]

taxas.previstas.ar <- Ar.Previsao(taxas.juro[intervalo.total, maturidade], 
                                  intervalo.passado,intervalo.futuro)

Multiplot.ts(data.frame(taxas.previstas.diebold[, maturidade],valores.reais[, maturidade],taxas.previstas.corte, taxas.previstas.rw))
(Eqm.diebold <- mean((taxas.previstas.diebold[, maturidade] - valores.reais[, maturidade])^2))
(Eqm.corte <- mean((taxas.previstas.corte - valores.reais[, maturidade])^2))
(Eqm.rw <- mean((taxas.previstas.rw - valores.reais[, maturidade])^2))

title(intervalo.total)




i <- 200
retirar <- taxas.juro[, 1]
curvas <- PreparaCurvasCorte(taxas.juro,percentual.testar,c(1,2500),s=125,maturidade=2,retirar=retirar)
valor.real <- PreparaValorFuturo(taxas.juro,percentual.testar,c(1,2500),s=61,maturidade=2)
# View(curvas$passado.learn)
# View(curvas$passado.test)
Distancias <- SemimetricPCA(curvas$passado.learn,curvas$passado.test,5)
# Distancias <- SemimetricDeriv(curvas$passado.learn,curvas$passado.test,0,5,c(0,10))
#View(Distancias)
# Multiplot.ts(t(curvas$passado.learn))
# Multiplot.ts(t(curvas$passado.test))
dist <- 1

# i <- 400
# plot(1:20,taxas.juro[i, ])
# plot(tempo.maturidades,taxas.juro[i, ])
tmp <- as.data.frame(cbind(Distancias=Distancias[which(Distancias[, i] < dist), i],
             Resposta = curvas$futuro.learn[which(Distancias[, i] < dist)],
             Diff.resp = abs(curvas$futuro.learn[which(Distancias[, i] < dist)] - 
                               curvas$futuro.learn[i]),
             Indice = which(Distancias[, i] < dist)))

Resposta <- curvas$futuro.learn + curvas$retirar.learn  # Resposta x Distancia
# Resposta <-  curvas$retirar.learn # Random walk

# compara a relação entre as respostas dos testes com a das respostas
# tmp2 <- data.frame(Distancias = Distancias[, i],
#                             Resp = Resposta,
#                             Diff.resp = abs(Resposta - valor.real[i]))
xyplot( Diff.resp ~ Distancias , data = tmp)
xyplot( Diff.resp ~ Distancias | equal.count(tmp2$Indice, number = 10), data = tmp2)

summary(lm(Diff.resp ~ Distancias,tmp2))

#View(tmp2)

#View(data.frame(Distancias=Distancias[, i],
#           Resp = Resposta,
#           Diff.resp = abs(Resposta - valor.real[i])))



#hist(tmp2$Diff.resp)
# 
# View(tmp[order(tmp$Distancias), ])
# data.plotar <- t(curvas$passado.learn[which(as.logical((Distancias[, i] < dist) *   
#                                                          (tmp2$Diff.resp > 2))), ])
# plot.ts(curvas$passado.test[i, ],ylim=c(min(data.plotar), 
#                                         max(data.plotar)))
#Multiplot.ts(data.plotar)
# 