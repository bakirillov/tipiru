library(psy)
library(plyr)
library(boot)
library(psych)
library(dplyr)
library(lavaan)
library(corrplot)
library(PerformanceAnalytics)

score <- function(v) {
  leikert.table <- data.frame(5,4,3,2,1)
  colnames(leikert.table) <- c(-2,-1,0,1,2)
  return(leikert.table[[as.character(v)]])
}

preprocess <- function(d) {
  for (a in colnames(d)[grep("X[0-9]+", colnames(d))]) {
    d[[a]] <- d[[a]] - 3
    d[[a]] <- sapply(d[[a]], score)
  }
  for (a in 2:11) {
    colnames(d)[a] <- paste("TIPI", toString(a-1), sep="")
  }
  colnames(d)[colnames(d) == "Ваш.пол."] <- "Sex"
  colnames(d)[colnames(d) == "Ваш.возраст."] <- "Age"
  d$Sex <- as.character(d$Sex)
  d$Sex[d$Sex == "Мужской"] <- "M"
  d$Sex[d$Sex == "Женский"] <- "F"
  d$Sex <- as.factor(d$Sex)
  return(d)
}

count5PFQ <- function(data) {
  r <- data
  r$b5.act.pass <- data$X1 + data$X6 + data$X11
  r$b5.dom.sub <- data$X16 + data$X21 + data$X26
  r$b5.jov.recl <- data$X31 + data$X36 + data$X41
  r$b5.seek.avoid <- data$X46 + data$X51 + data$X56
  r$b5.attr.avoid <- data$X61 + data$X66 + data$X71
  r$b5.warmth.indf <- data$X7 + data$X12 + data$X2
  r$b5.coop.conc <- data$X17 + data$X22 + data$X27
  r$b5.trst.susp <- data$X32 + data$X37 + data$X42
  r$b5.undst.not <- data$X47 + data$X52 + data$X57
  r$b5.resp.oth.slf <- data$X62 + data$X67 + data$X72
  r$b5.neat.not <- data$X3 + data$X8 + data$X13
  r$b5.ins.ww <- data$X18 + data$X23 + data$X28
  r$b5.resp.irr <- data$X33 + data$X38 + data$X43
  r$b5.slfcnt.imp <- data$X48 + data$X53 + data$X58
  r$b5.4sght.cfr <- data$X63 + data$X68 + data$X73
  r$b5.anx.cfr <- data$X4 + data$X9 + data$X14
  r$b5.tns.rlx <- data$X19 + data$X24 + data$X29
  r$b5.dpr.comf <- data$X34 + data$X39 + data$X44
  r$b5.scr.ssf <- data$X49 + data$X54 + data$X59
  r$b5.lab.stab <- data$X64 + data$X69 + data$X74
  r$b5.cur.cons <- data$X5 + data$X10 + data$X15
  r$b5.drm.real <- data$X20 + data$X25 + data$X30
  r$b5.art.not <- data$X35 + data$X40 + data$X45
  r$b5.sens.not <- data$X50 + data$X55 + data$X60
  r$b5.pls.rig <- data$X65 + data$X70 + data$X75
  r$E.5PFQ <- r$b5.dom.sub + r$b5.act.pass + r$b5.jov.recl + r$b5.seek.avoid + r$b5.attr.avoid
  r$A.5PFQ <- r$b5.warmth.indf + r$b5.coop.conc + r$b5.trst.susp + r$b5.undst.not + r$b5.resp.oth.slf
  r$C.5PFQ <- r$b5.neat.not + r$b5.ins.ww + r$b5.resp.irr + r$b5.slfcnt.imp + r$b5.4sght.cfr
  r$ES.5PFQ <- r$b5.anx.cfr + r$b5.tns.rlx + r$b5.dpr.comf + r$b5.scr.ssf + r$b5.lab.stab
  r$O.5PFQ <- r$b5.cur.cons + r$b5.pls.rig + r$b5.sens.not + r$b5.art.not + r$b5.drm.real
  return(r)
}

countTIPI <- function(data) {
  upsideDown <- function(x) {rev(list(1,2,3,4,5,6,7))[[x]]}
  r <- data.frame(
    X=data$X., TIPI1=data$TIPI1, TIPI2=data$TIPI2, TIPI3=data$TIPI3,
    TIPI4=data$TIPI4, TIPI5=data$TIPI5, TIPI6=data$TIPI6, TIPI7=data$TIPI7,
    TIPI8=data$TIPI8, TIPI9=data$TIPI9, TIPI10=data$TIPI10
  )
  r$E.TIPI <- 0.5*(r$TIPI1 + sapply(r$TIPI6, upsideDown))
  r$A.TIPI <- 0.5*(r$TIPI7 + sapply(r$TIPI2, upsideDown))
  r$C.TIPI <- 0.5*(r$TIPI3 + sapply(r$TIPI8, upsideDown))
  r$ES.TIPI <- 0.5*(r$TIPI9 + sapply(r$TIPI4, upsideDown))
  r$O.TIPI <- 0.5*(r$TIPI5 + sapply(r$TIPI10, upsideDown))
  return(r)
}

getOmega2 <- function(first, second) {
  d <- cbind(scale(first), 1-scale(second))
  return(omega(d, nfactors = 1)[["omega.tot"]])
}

getAlpha <- function(first, second) {
  d <- cbind(scale(first), 1-scale(second))
  return(cronbach(d)[["alpha"]])
}

toCrb <- function (first, second) {
  return(cbind(scale(first), 1-scale(second)))
}

bootAlpha <- function(d) {
  cronbach.boot <- function (data, x) {cronbach(data[x,])[[3]]}
  res <- boot(d, cronbach.boot, 1000)
  q <- quantile(res$t, c(0.025, 0.975))
  bca <- boot.ci(res, type="bca")[["bca"]]
  return(rbind(q[["2.5%"]], q[["97.5%"]], bca[[4]], bca[[5]]))
}
