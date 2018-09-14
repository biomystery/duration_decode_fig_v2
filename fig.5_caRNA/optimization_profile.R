require(optimization)
hi <- function(x){(x[1]**2 + x[2] - 11)**2 + (x[1] + x[2]**2 -7)**2}
optim_nm(fun = hi, k = 2) # a direct search algorithm 
optim_sa(fun = hi, start = c(runif(2, min = -1, max = 1)),
         trace = FALSE,
         lower = c(-4, -4),
         upper = c(4, 4),
         control = list(dyn_rf = FALSE,
                        rf = 1.2,
                        t0 = 10,
                        nlimit = 100,
                        r = 0.6,
                        t_min = 0.1
         )
)


##### Himmelblau's function
# minimum at f(3,2) = 0
# f(-2.805, -3.1313) = 0
# f(-3.779, -3.283) = 0
#f(3.5844, -1.848) = 0
H <- function(x){
  (x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
}
Output <- optim_nm(H, k = 2, trace = TRUE)
plot(Output) 
plot(Output, 'contour',lower = c(-4, -4), upper = c(4, 4))

out_sa <- optim_sa(fun = H, start = c(runif(2, min = -1, max = 1)),
                   trace = TRUE, lower = c(-4, -4) ,upper=c(4, 4),
                   control = list(t0 = 1000, nlimit = 1500,r = 0.8))
plot(out_sa) # simulated annealing
plot(out_sa, type = "contour")


# nlminb method -----------------------------------------------------------
# use capture.output to capture traces 
x <- rnbinom(100, mu = 10, size = 10)
hdev <- function(par)
  -sum(dnbinom(x, mu = par[1], size = par[2], log = TRUE))
output <-capture.output(nlminb(c(9, 12), hdev,control=list(trace=T)),type = 'output',split = T)
par(mfrow=c(3,1))
sapply(1:3, function(i) plot(Output$trace[,c(1,i+1)],type = 'b'))
par(mfrow=c(1,1))


# real time plotting ------------------------------------------------------
n=1000
df=data.frame(time=1:n,y=runif(n))
window=100
for(i in 1:(n-window)) {
  flush.console()
  #plot(df$time,df$y,type='l',xlim=c(i,i+window))
  plot(df$time,df$y,type='l',xlim=c(0,1000))
  Sys.sleep(.09)
}



# annimation with ggplot --------------------------------------------------
library(animation)
library(ggplot2)
# your data
n  <- 1000
df <- data.frame(time=1:n,y=runif(n))
window <- 100
# create the animation
saveHTML({
  for(i in 1:(n-window)) {
    print(ggplot(df) + geom_line(aes(x=time, y=y), size=0.7) + xlim(i,i+window))
  }
})


# gif production ----------------------------------------------------------
ani.options(interval = 0.05)
saveGIF(runa())
