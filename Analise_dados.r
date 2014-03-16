

# Resultados gerais -------------------------------------------------------

tabela <- TabelaEQM(simulacoes,vetor.maturidade,horizonte=12)
heatmap(tabela,scale="row", Colv=NA, Rowv=NA)
View(tabela)
PlotarSeries(simulacoes,maturidade=1,horizonte=3,subconjunto=c(1,20))



tabelao.agrupado <- TabelaoAgrupado(simulacoes,vetor.maturidade,horizonte=5,vetor.intervalos)
View(tabelao.agrupado); rowSums(tabelao.agrupado)
tabelao.agrupado <- tabelao.agrupado/tabelao.agrupado[, 3]
heatmap(tabelao.agrupado[c(-8,-9), ], scale="column", Colv=NA, Rowv=NA) # tira as estimacoes do que não retiram nada

tabelao.agrupado <- TabelaoAgrupado(simulacoes,vetor.maturidade,horizonte=25,vetor.intervalos)
View(tabelao.agrupado); rowSums(tabelao.agrupado)
heatmap(tabelao.agrupado[c(-8,-9), ], scale="column", Colv=NA, Rowv=NA) # tira as estimacoes do que não retiram nada

tabelao.agrupado <- TabelaoAgrupado(simulacoes,vetor.maturidade,horizonte=50,vetor.intervalos)
View(tabelao.agrupado); rowSums(tabelao.agrupado)
heatmap(tabelao.agrupado[c(-8,-9), ], scale="column", Colv=NA, Rowv=NA) # tira as estimacoes do que não retiram nada

GiacominiWhite(simulacoes,benchmarking=3,valor.real=12,tamanho=800,horizonte=25,vetor.maturidade,vetor.intervalos)


# Análise detalhada ----------------------------------------------------


# investiga qualquer item (heatmap e tabela)
tabelao <- Tabelao(simulacoes,vetor.maturidade,vetor.horizontes=c(5),vetor.intervalos)
tabelao <- tabelao[, c(-8,-9)] # tira as estimacoes do que não retiram nada, pois o resultado é muito ruim
require(stringr)
par.ordenar <- 1; tabelao <- tabelao[order(str_split_fixed(rownames(tabelao), "[|]",3)[, par.ordenar]), ] # par.ordenar define qual por qual dos 3 campos que será ordenada as estimações
tabelao <- CombinarPrevisoes(tabelao,"crt.p5.rmat6","rw")
# tabelao <- CombinarPrevisoes(tabelao,"corte.deriv.1","rw")
heatmap(tabelao, scale="row",Rowv=NA) # scale="row",Rowv=NA
require(lattice)
levelplot(tabelao, scale=list(x=list(rot=45)))
summary(tabelao)
View(tabelao)


