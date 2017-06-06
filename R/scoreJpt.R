scoreJpt <- function(dat,id,listUp,listDown){
  allsc <- allpv <- allgp <- NULL

  for(i in 1:ncol(dat)){

    Fn <- ecdf(dat[,i])

    ecdfup <- Fn(dat[,i])[id %in% listUp]
    ecdfdn <- Fn(dat[,i])[id %in% listDown]

    restest <- ks.test(x=ecdfup, y=ecdfdn)
    score <- restest$statistic
    pv <- restest$p.value
    gp <- NA

    if(mean(ecdfup) < mean(ecdfdn)) {
      score <- -score
      if(pv < 0.005){ gp <- "Epi" }else{gp <- "Epi-Int"}
    }else if (mean(ecdfup) > mean(ecdfdn)){
      if(pv < 0.005){ gp <- "Mes" }else{gp <- "Mes-Int"}
    }

    allsc <- c(allsc,score)
    allpv <- c(allpv,pv)
    allgp <- c(allgp,gp)
  }

  res <- cbind(colnames(dat),allsc,allpv,allgp)
  return(res)
}

