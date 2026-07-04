#par.old=FLPar(Linf=238.59,K=0.185,t0=-1.404,a=1.76e-3,b=3.378,l50=142,
#          m1=20/c(mean(stock.wt(swom)["3"])),m2=FLPar(-0.309))
library(FLCore)
library(FLife)

swom <- readRDS("Robj/swo_stk.rds")
par=FLPar(Linf=238.59,K=0.185,t0=-1.404,a=9.62e-3,b=3.06,l50=132,
          m1=20/c(mean(stock.wt(swom)["3"])),m2=FLPar(-0.309))

m=FLife:::lorenzen(stock.wt(swom),par[c("m1","m2")])
par["m1"]=0.2*par["m1"]/mean(m["3"])  

m(swom)=FLife:::lorenzen(stock.wt(swom),par[c("m1","m2")])

lor_m=as.data.frame(swom@m)
ggplot(lor_m,(aes(age,data)))+
  geom_line(stat="identity")

ggplot(lor_m,(aes(age,data)))+
  geom_smooth(method="loess")+
  geom_point()

saveRDS(swom, file = "Robj/swo_stk_loren.rds")
catch.n(swom)[catch.n(swom) == 0]
catch.n(swom)
