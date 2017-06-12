library(spatstat)
library(rstan)
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
source("utils.r")

load("WashingtonDC-agls.rdata")

args <- commandArgs(trailingOnly = TRUE)

print(args)
start_date = as.Date(args[1]) 
end_date = as.Date(args[2]) 
bw_space = as.numeric(args[3])
bw_time = as.numeric(args[4])
separable = as.numeric(args[5]) == 1
remove.duplicates = as.numeric(args[6]) == 1

## Preprocess the data based on the command line arguments
keep = which(data$date >= start_date & data$date < end_date & data$holiday == F)
xyt = xyt[keep,]
xy= xy[keep,]
nrow(xyt)

if(remove.duplicates) {
  time_dist = dist(xyt[,3])
  space_dist = dist(xyt[,1:2])
  library(igraph)
  ig = graph.adjacency(as.matrix(space_dist) < .1 & as.matrix(time_dist) <= .5  & lower.tri(space_dist))
  n.orig = nrow(xyt)
  clust = clusters(ig)$membership
  keep = !duplicated(clust)
  xyt = xyt[keep,]
  xy = xy[keep,]
}

xyt$T = xyt$T / 24
xy.area = area.owin(my.owin)
time.window = as.numeric(end_date-start_date)+1

xyt$hour = round(xyt$T*24) %% 24 + 1
hourly.density = table(factor(xyt$hour,levels=1:24)) / nrow(xyt)
phat.hour = hourly.density[xyt$hour]

phat.space = density.ppp(xy, at="points", sigma=bw_space) / nrow(xyt)

if(separable) {
	muhat = intensity.estimate.separable(xyt,bw_space,bw_time)
} else {
	muhat = intensity.estimate(xyt,bw_space,bw_time)
}
muSTintegral = nrow(xyt)

## Inference
m = stan_model("hawkes-model.stan")
data = list(Space=xyt[,1:2],Time=xyt[,3]-min(xyt[,3]), muST=muhat, muSTintegral=muSTintegral, 
	    n=nrow(xyt),time_window=time.window, space_window=xy.area, hour=round(xyt$T*24) %% 24+1)

fit = sampling(m,data,warmup=100,iter=200,chains=4) 
print(fit,c("lengthscaleS","lengthscaleT","a","mu"))
print("Time decay:")
out = extract(fit)
print(quantile(1/out$lengthscaleT,c(.025,.1,.5,.9,.975)))

fname = sprintf("data/hawkes-model-%s-%s-%.02f-%.02f-%d-%d.rdata", start_date,end_date,bw_space,bw_time,separable,remove.duplicates)
save(l,fit,out,data,file=fname)
print(fname)
print("Background attribution: "); quantile(apply(out$background,1,mean),c(.025,.5,.975))
