(
  WD <- getwd()
)
if (!is.null(WD)) setwd(WD)

dfAS <- read.csv("WSE-ASSECOPOL.csv",sep = ",")
dfAS <- dfAS[nrow(dfAS):1,]
dfCD <- read.csv("WSE-CDPROJEKT.csv",sep = ",")
dfCD <- dfCD[nrow(dfCD):1,]
dfCY <- read.csv("WSE-CYFRPLSAT.csv",sep = ",")
dfCY <- dfCY[nrow(dfCY):1,]

stopaPR <- function(.df,.df1,.df2){
  VolatilityAS <- NULL
  VolatilityCD <- NULL
  VolatilityCY <- NULL
  for(i in 1:(nrow(.df)-1)){
    VolatilityAS[i] <- ((.df$Close[i+1]-.df$Close[i])/(.df$Close[i]))
    VolatilityCD[i] <- ((.df1$Close[i+1]-.df1$Close[i])/(.df1$Close[i]))
    VolatilityCY[i] <- ((.df2$Close[i+1]-.df2$Close[i])/(.df2$Close[i]))
  }
  #Volatility[nrow(.df)] <- NA
  final <- cbind(VolatilityAS,VolatilityCD,VolatilityCY)
  return(final)
}

EWMA <- function(.df){
  .varAS <- NULL
  .varCD <- NULL
  .varCY <- NULL
  .covAS.CD <- NULL
  .covAS.CY <- NULL
  .covCD.CY <- NULL
  lambda <- 0.94
  .varAS[1] <- var(.df$VolatilityAS[1:250])
  .varCD[1] <- var(.df$VolatilityCD[1:250])
  .varCY[1] <- var(.df$VolatilityCY[1:250])
  .covAS.CD[1] <- cov(.df$VolatilityAS[1:250],.df$VolatilityCD[1:250])
  .covAS.CY[1] <- cov(.df$VolatilityAS[1:250],.df$VolatilityCY[1:250])
  .covCD.CY[1] <- cov(.df$VolatilityCD[1:250],.df$VolatilityCY[1:250])
  for(i in 1:nrow(.df)){
    .varAS[i+1] <- .varAS[i]*lambda + ((1-lambda)*.df$VolatilityAS[i]^2)
    .varCD[i+1] <- .varCD[i]*lambda + ((1-lambda)*.df$VolatilityCD[i]^2)
    .varCY[i+1] <- .varCY[i]*lambda + ((1-lambda)*.df$VolatilityCY[i]^2)
    .covAS.CD[i+1] <- .covAS.CD[i]*lambda + ((1-lambda)*.df$VolatilityAS[i]*.df$VolatilityCD[i])
    .covAS.CY[i+1] <- .covAS.CY[i]*lambda + ((1-lambda)*.df$VolatilityAS[i]*.df$VolatilityCY[i])
    .covCD.CY[i+1] <- .covCD.CY[i]*lambda + ((1-lambda)*.df$VolatilityCD[i]*.df$VolatilityCY[i])
  }
  final <- cbind(.varAS,.varCD,.varCY,.covAS.CD,.covAS.CY,.covCD.CY)
  final <- as.data.frame(final)
  return(final)
}

mVaR <- function(.df){
  VaR99 <- NULL
  ES99 <- NULL
  VaR95 <- NULL
  ES95 <- NULL
  .corAS.CD <- NULL
  .corAS.CY <- NULL
  .corCD.CY <- NULL
  .row1 <- c(0.33,0.33,0.33)
  for(i in 2:nrow(.df)){
    .row2 <- c(.df$.varAS[i],.df$.covAS.CD[i],.df$.covAS.CY[i])
    .row3 <- c(.df$.covAS.CD[i],.df$.varCD[i],.df$.covCD.CY[i])
    .row4 <- c(.df$.covAS.CY[i],.df$.covCD.CY[i],.df$.varCY[i])
    macierz <- rbind(.row2,.row3,.row4)
    macierz1 <- .row1*macierz
    temp <- c(sum(macierz1[,1]),sum(macierz1[,2]),sum(macierz1[,3]))
    war.portfela <- sum(.row1*temp)
    SD.portfela <- sqrt(war.portfela)
    VaR99[i-1] <- 2.33*SD.portfela*100
    VaR95[i-1] <- 1.64*SD.portfela*100
    ES99[i-1] <- SD.portfela*(exp(-(2.33^2)/2))/(sqrt(2*pi)*0.01)
    ES95[i-1] <- SD.portfela*(exp(-(1.64^2)/2))/(sqrt(2*pi)*0.05)
    .SD <- c(sqrt(macierz[1,1]),sqrt(macierz[2,2]),sqrt(macierz[3,3]))
    .corAS.CD[i-1] <- macierz[1,2]/(.SD[1]*.SD[2])
    .corAS.CY[i-1] <- macierz[1,3]/(.SD[1]*.SD[3])
    .corCD.CY[i-1] <- macierz[2,3]/(.SD[2]*.SD[3])
  }
  final <- cbind(VaR99,ES99,VaR95,ES95,.corAS.CD,.corAS.CY,.corCD.CY)
  final <- as.data.frame(final)
  return(final)
}
test <- stopaPR(dfAS,dfCD,dfCY)
test <- as.data.frame(test)
zxc <- EWMA(test)
f <- mVaR(zxc)

plot(f$VaR99, type = "l")
lines(f$VaR95,type = "l", col = "red")

plot(f$ES99, type = "l")
lines(f$ES95,type = "l", col = "red")

plot(f$.corAS.CD, type = "l")
lines(f$.corAS.CY,type = "l", col = "red")
lines(f$.corCD.CY,type = "l", col = "green")