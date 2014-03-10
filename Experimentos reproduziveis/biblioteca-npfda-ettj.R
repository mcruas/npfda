######################### Trabalho teve como ponto de partida como ponto de 
#partida os arquivos disponíveis em 
#http://www.math.univ-toulouse.fr/staph/npfda/
#
# Um pequeno tutorial sobre como utilizar as funções abaixo pode ser obtido no
# endereço abaixo:
# http://www.rpubs.com/mcruas/npfda-ettj
#########################

library("stringr")

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

PreparaCurvasCorte  <- function(base,maturidade,intervalo.passado,intervalo.futuro,todas.maturidades = NULL,retirar = NULL) {
  ########################################################################
  # Essa função retorna uma variável com as curvas necessárias para fazer a estimação.
  # 
  # ARGS
  #   base: uma base de dados; cada linha representa uma série temporal
  #   percentual.testar: o percentual da amostra que deve ser utilizado para testar 
  #         a capacidade de previsão do modelo
  #   intervalo.passado: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo do treino
  #   intervalo.futuro: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo a ser previsto
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
  horizonte <- intervalo.futuro[1] - intervalo.passado[2]
  learning <- intervalo.passado - c(0,horizonte)
  testing <- intervalo.futuro - horizonte
  curvas <- NULL
  curvas$retirar.learn <- retirar[learning[1]:learning[2]]
  curvas$retirar.test <- retirar[testing[1]:testing[2]]
  curvas$passado.learn <- as.matrix(base[learning[1]:learning[2], ] - retirar[learning[1]:learning[2]])
  curvas$passado.test <- as.matrix(base[testing[1]:testing[2], ] - retirar[testing[1]:testing[2]])
  curvas$futuro.learn <- as.matrix(base[learning[1]:learning[2]+horizonte,maturidade] - retirar[learning[1]:learning[2]])
  curvas$estimacao  <- NULL
  curvas$tipo  <- "corte"
  class(curvas) <- "fdaCorte"
  return(curvas)
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
