###On settlement date Dec 13, 2021, the following interest rate swap was initiated. 
###Fixed rate:  0.853% 
###Floating rate:  3 month LIBOR 
###Notional principal: $100 million 
###Maturity:  2 years 
###Frequency:  quarterly payments 
###Day count:  actual/360 
###Settlement:  2021-12-15 

###The first floating rate payment is based on the 3-month LIBOR rate of 0.20275%  

###We want to value the swap based on the following information: 
###LIBOR rates, for settlement on 2/2/2022 (libor.csv) 
###ED futures prices for settlement on 2/2/2022 (edfut.csv) 

### start

###input the files
libor<-read.csv("libor.csv")
edf<-read.csv("edfut.csv")
libor
edf
###functions
DATE <- function(yyyy,mm,dd) {
  dte  <- as.Date(sprintf("%i-%i-%i",yyyy,mm,dd),format="%Y-%m-%d")
  return(dte)
}
as.Date2 <- function(x) {
  tryfmt <- c("%Y-%m-%d","%m/%d/%Y","%Y/%m/%d","%b %d,%Y")
  return(as.Date(x,tryFormats=tryfmt))
}
###create data.table for swap payment dates
swap <- data.table(period=c(0:8),
                   pay.date=seq(DATE(2021,12,15),DATE(2023,12,15),by="3 month"))
swap[, ndays := c(0,as.numeric(diff(swap$pay.date))) ]
knitr::kable(swap)
###interpolate libor rate
settle <- DATE(2021,12,15)
names(libor) <- c("term","spot")
swap[, L3 := NA]
swap$L3[1] <- 0.0020275
knitr::kable(swap)

names(edf) <- c("Expmon","FutMatDt","FutPr")
nedf <- nrow(edf)

# remove futures that do not mature in Mar/Jun/Sep/Dec
edf$Expmon<-as.character(edf$Expmon)
edf <- edf[substr(edf$Expmon,5,6) %in% c("03","06","09","12"),]
edf$FutMatDt <- as.Date2(edf$FutMatDt)
edf$start.date <- as.Date2(edf$FutMatDt)+2          # T+2 settlement
edf$forw <- 1-edf$FutPr/100 

# include today's 3 month libor (spot) rate
tmp <- data.table(edf[,c(1,4,5)])
setorderv(tmp,c("start.date"))
tmp
# interpolate 3 month forward rates
L3 <- as.data.frame(spline(x=tmp$start.date,y=tmp$forw,
                           xout=swap$pay.date,method="natural"))$y
swap$L3[2:nrow(swap)] <- L3[2:length(L3)]
knitr::kable(swap)
###find cash flows
notional   <-100
swap$cf.float <- 0
swap$cf.fixed <- 0
swap.rate.01 <- 0.00853
for (i in 2:nrow(swap)) {
  swap$cf.float[i] <- swap$ndays[i]*swap$L3[i-1]/360 * notional
  swap$cf.fixed[i] <- swap$ndays[i]*swap.rate.01/360 * notional
}
knitr::kable(swap)
###interpolate spot rate from 2022-02-02 to 2022-03-15
libor<-cbind(libor,start.date = 
               c("2022-02-03","2022-03-02","2022-05-02",
                 "2022-08-02","2023-02-02"))

libor$start.date <-as.Date(libor$start.date)

spot<-spline(x=libor$start.date,y=libor$spot,xout=DATE(2022,03,15),
             method="natural")$y
swap$L3[1]<-spot

###find discount factor
swap$disfac <- 1
for (i in 2:nrow(swap)) {
  nday <- as.numeric(swap$pay.date[i]-swap$pay.date[i-1])
  swap$disfac[i] <- swap$disfac[i-1] / (1+swap$L3[i-1]*nday/360)
}

knitr::kable(swap)
###present value
PV_float <- sum(swap$disfac*swap$cf.float)
PV_fixed <- sum(swap$disfac*swap$cf.fixed)
PV_fixed  - PV_float 



### end
